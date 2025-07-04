#!/usr/bin/env python3
"""Calculate frequency and time offsets between two parties.

The timestamps can either be as 'a1' timestamp files (per readevents),
or as 'T1'/'T2' epoch files (per qcrypto).

Since the time offset will evolve with clock skew, the reference time is defined
as the common starting time of both parties, and the time offset at said reference.

The results are written to stdout, in the order -
df (absolute), \\t, dt (ns) - where 'df' is the compensating
frequency offset and 'dt' the compensating timing offset for the target
timestamps, i.e. '(TARGET - dt) / (1 + df)'.

To scan for possible precompensations, use the '--precomp-enable' flag.
For more detailed options, use multiple help flags.

Changelog:
    2023-01-09, Justin: Refactoring from fpfind.py
    2023-01-31, Justin: Formalize interface for fpfind.py
"""

import functools
import inspect
import os
import sys
import time
from pathlib import Path

import configargparse
import numpy as np

from fpfind.lib._logging import get_logger, set_logfile, verbosity2level
from fpfind.lib.constants import (
    EPOCH_LENGTH,
    MAX_FCORR,
    NTP_MAXDELAY_NS,
    PeakFindingFailed,
)
from fpfind.lib.parse_timestamps import read_a1, read_a1_start_end
from fpfind.lib.utils import (
    ArgparseCustomFormatter,
    generate_fft,
    get_first_overlapping_epoch,
    get_statistics,
    get_timestamp_pattern,
    get_timing_delay_fft,
    get_xcorr,
    normalize_timestamps,
    parse_docstring_description,
    round,
    slice_timestamps,
)

logger = get_logger(__name__, human_readable=True)

# Allows quick prototyping by interrupting execution right before FFT
# Trigger by running 'python3 -i -m fpfind.fpfind --config ... --experiment'.
# For internal use only.
_ENABLE_BREAKPOINT = False

# Disables resolution doubling during initial peak search, used to converge
# on specific fpfind parameters for peak search. This procedure is a relatively
# cheap measure to increase peak finding yield, and thus is unlikely needed
# in production.
# For internal use only.
_DISABLE_DOUBLING = False

# Allow fpfind to automatically perform recovery actions when
# insufficient coincidences is assumed to be in the cross-correlation calculation
# For internal use only.
_NUM_WRAPS_LIMIT = 0

# Controls learning rate, i.e. how much to decrease resolution by
RES_REFINE_FACTOR = np.sqrt(2)

# Set lower limit to frequency compensation before disabling
TARGET_DF = 1e-10

PROFILE_LEVEL = 0


def profile(f):
    """Performs a simple timing profile.

    Modifies the logging facility to log the name and lineno of the function
    being called, instead of the usual encapsulating function.

    Possibly thread-safe, but may not yield the correct indentation levels
    during execution, in that situation.
    """

    # Allows actual function to be logged rather than wrapped
    caller = inspect.getframeinfo(inspect.stack()[1][0])
    extras = {
        "_funcname": f"[{f.__name__}]",
        "_filename": os.path.basename(caller.filename),
        # "_lineno": caller.lineno,  # this returns the @profile lineno instead
    }

    @functools.wraps(f)
    def wrapper(*args, **kwargs):
        global PROFILE_LEVEL
        pad = "  " * PROFILE_LEVEL
        logger.debug(
            "%sSTART PROFILING",
            pad,
            stacklevel=2,
            extra=extras,
        )
        PROFILE_LEVEL += 1

        try:
            start = time.time()
            result = f(*args, **kwargs)
        finally:
            end = time.time()
            elapsed = end - start
            logger.debug(
                "%sEND PROFILING: %.0f s",
                pad,
                elapsed,
                stacklevel=2,
                extra=extras,
            )
            PROFILE_LEVEL -= 1

        return result

    return wrapper


# Main algorithm
@profile
def time_freq(
    ats: list,
    bts: list,
    num_wraps: int,
    num_bins: int,
    resolution: float,
    target_resolution: float,
    threshold: float,
    separation_duration: float,
    quick: bool,
    do_frequency_compensation: bool = True,
):
    """Perform the actual frequency compensation routine.

    Timestamps must already by normalized to starting time of 0ns. Whether
    this starting time is also common reference between 'ats' and 'bts' is up to
    implementation.

    Args:
        ats: Timestamps of reference side, in units of ns.
        bts: Timestamps of compensating side, in units of ns.
        num_wraps: Number of cross-correlations to overlay, usually 1.
        num_bins: Numbers of bins to use in the FFT.
        resolution: Initial resolution of cross-correlation.
        target_resolution: Target resolution desired from routine.
        threshold: Height of peak to discriminate as signal, in units of dev.
        separation_duration: Separation of cross-correlations, in units of ns.
    """
    end_time = min(ats[-1], bts[-1])
    duration = num_wraps * resolution * num_bins
    logger.debug("  Performing peak searching...")
    logger.debug(
        "    Parameters:",
        extra={
            "details": [
                f"Bins: 2^{np.int32(np.log2(num_bins)):d}",
                f"Bin width: {resolution:.0f}ns",
                f"Number of wraps: {num_wraps:d}",
                f"Target resolution: {target_resolution}ns",
            ]
        },
    )

    # Refinement loop, note resolution/duration will change during loop
    dt = 0
    f = 1
    curr_iteration = 1

    while True:
        # Dynamically adjust 'num_wraps' based on current 'resolution',
        # avoids event overflow/underflow.
        #
        # Previous implementation tried to use as many events as possible
        # This became a problem when fine-tuning the timing resolution, due
        # to the presence of multiple peaks under a clock skew, which will now
        # be resolvable. Instead, it is sufficient to force the same 'num_wraps'.
        #
        # TODO(2024-02-14): Account for the separation time when calculating
        #                   'max_wraps'.
        max_wraps = np.floor(end_time / (resolution * num_bins))
        _num_wraps = min(num_wraps, max(max_wraps, 1))  # bounded by 1 <= k <= k_max
        _duration = _num_wraps * resolution * num_bins
        logger.debug(
            "    Iteration %s (r=%.1fns, k=%d)",
            curr_iteration,
            resolution,
            _num_wraps,
        )

        # Perform cross-correlation
        logger.debug(
            "      Performing earlier xcorr (range: [0.00, %.2f]s)",
            _duration * 1e-9,
        )
        ats_early = slice_timestamps(ats, 0, _duration)
        bts_early = slice_timestamps(bts, 0, _duration)
        afft = generate_fft(ats_early, num_bins, resolution)
        bfft = generate_fft(bts_early, num_bins, resolution)
        ys = get_xcorr(afft, bfft)

        # Calculate timing delay
        _dtype = np.int32  # signed number needed for negative delays
        if num_bins * resolution > 2147483647:  # int32 max
            _dtype = np.int64
        xs = np.arange(num_bins, dtype=_dtype) * resolution
        dt1 = get_timing_delay_fft(ys, xs)[0]  # get smaller candidate
        sig = get_statistics(ys, resolution).significance

        # Confirm resolution on first run
        # If peak finding fails, 'dt' is not returned here:
        # guaranteed zero during first iteration
        while curr_iteration == 1:
            logger.debug(
                "        Peak: S = %.3f, dt = %sns (resolution = %.0fns)",
                sig,
                dt1,
                resolution,
            )

            # Deviation zero, due flat cross-correlation
            # Hint: Increase number of bins
            if sig == 0:
                raise PeakFindingFailed(
                    "Bin saturation  ",
                    significance=sig,
                    resolution=resolution,
                    dt1=dt1,
                )

            # If peak rejected, merge contiguous bins to double
            # the resolution of the peak search
            if sig >= threshold:
                logger.debug("          Accepted")
                break
            else:
                logger.debug("          Rejected")
                if _DISABLE_DOUBLING:
                    raise PeakFindingFailed(
                        "Low significance",
                        significance=sig,
                        resolution=resolution,
                        dt1=dt1,
                    )
                ys = np.sum(ys.reshape(-1, 2), axis=1)
                resolution *= 2

            # Catch runaway resolution doubling, limited by
            if resolution > 1e4:
                raise PeakFindingFailed(
                    "Resolution OOB  ",
                    significance=sig,
                    resolution=resolution,
                    dt1=dt1,
                )

            # Recalculate timing delay since resolution changed
            xs = np.arange(len(ys)) * resolution
            dt1 = get_timing_delay_fft(ys, xs)[0]
            sig = get_statistics(ys, resolution).significance

        # Calculate some thresholds to catch when peak was likely not found
        buffer = 1
        threshold_dt = buffer * max(abs(dt), 1)
        threshold_df = buffer * max(abs(f - 1), 1e-9)

        # If the duration has too few coincidences, the peak may not
        # show up at all. A peak is likely to have already been found, if
        # the current iteration is more than 1. Attempt to retry with
        # larger 'num_wraps' if enabled in settings.
        if curr_iteration != 1 and abs(dt1) > threshold_dt:
            logger.warning(
                "      Interrupted due spurious signal:",
                extra={
                    "details": [
                        f"early dt       = {dt1:10.0f} ns",
                        f"threshold dt   = {threshold_dt:10.0f} ns",
                    ]
                },
            )
            if num_wraps < _NUM_WRAPS_LIMIT:
                num_wraps += 1
                logger.debug(
                    "      Reattempting with k = %d due missing peak.",
                    num_wraps,
                )
                continue

        # Catch if timing delay exceeded
        if abs(dt1) > NTP_MAXDELAY_NS:
            raise PeakFindingFailed(
                "Time delay OOB  ",
                significance=sig,
                resolution=resolution,
                dt1=dt1,
            )

        _dt1 = dt1
        df1 = 0  # default values
        if do_frequency_compensation:
            # TODO(2024-01-31):
            #     If signal too low, spurious peaks may occur. Implement
            #     a mechanism to retry to obtain an alternative value.
            logger.debug(
                "      Performing later xcorr (range: [%.2f, %.2f]s)",
                separation_duration * 1e-9,
                (separation_duration + _duration) * 1e-9,
            )
            ats_late = slice_timestamps(ats, separation_duration, _duration)
            bts_late = slice_timestamps(bts, separation_duration, _duration)
            _afft = generate_fft(ats_late, num_bins, resolution)
            _bfft = generate_fft(bts_late, num_bins, resolution)
            _ys = get_xcorr(_afft, _bfft)

            # Calculate timing delay for late set of timestamps
            xs = np.arange(num_bins) * resolution
            _dt1 = get_timing_delay_fft(_ys, xs)[0]
            df1 = (_dt1 - dt1) / separation_duration

            # Some guard rails to make sure results make sense
            # Something went wrong with peak searching, to return intermediate
            # results which are likely near correct values.
            if curr_iteration != 1 and (
                abs(_dt1) > threshold_dt or abs(df1) > threshold_df
            ):
                logger.warning(
                    "      Interrupted due spurious signal:",
                    extra={
                        "details": [
                            f"early dt       = {dt1:10.0f} ns",
                            f"late dt        = {_dt1:10.0f} ns",
                            f"threshold dt   = {threshold_dt:10.0f} ns",
                            f"current df     = {df1*1e6:10.4f} ppm",
                            f"threshold df   = {threshold_df*1e6:10.4f} ppm",
                        ]
                    },
                )
                if num_wraps < _NUM_WRAPS_LIMIT:
                    num_wraps += 1
                    logger.debug(
                        "      Reattempting with k = %d due missing peak.",
                        num_wraps,
                    )
                    continue

                break  # terminate if recovery not enabled

            # Catch if timing delay exceeded
            if abs(_dt1) > NTP_MAXDELAY_NS:
                raise PeakFindingFailed(
                    "Time delay OOB  ",
                    significance=sig,
                    resolution=resolution,
                    dt1=_dt1,
                )

        # Apply recursive relations
        #   A quick proof (note dt -> t):
        #     iter 0 ->  (T - t0)/f0          = T/f0    - (t0/f0)
        #     iter 1 -> ((T - t0)/f0 - t1)/f1 = T/f0/f1 - (t0/f0/f1 + t1/f1)
        #   We want:
        #     iter n ->  (T - t)/f            = T/f     - t/f
        #   Thus:
        #     f = f0 * f1 * ... * fn
        #     t = f * (t0/f0/f1/.../fn + t1/f1/.../fn + ... + tn/fn)
        #       = t0 + t1*f0 + t2*f0*f1 + ... + tn*f0*f1*...*(fn-1)
        #   Recursive relation:
        #     f' = f * fn
        #     t' = f * tn + t,  i.e. use old value of f
        dt += f * dt1
        f *= 1 + df1
        logger.debug(
            "      Calculated timing delays:",
            extra={
                "details": [
                    f"early dt       = {dt1:10.0f} ns",
                    f"late dt        = {_dt1:10.0f} ns",
                    f"accumulated dt = {dt:10.0f} ns",
                    f"current df     = {df1*1e6:10.4f} ppm",
                    f"accumulated df = {(f-1)*1e6:10.4f} ppm",
                ]
            },
        )

        # Throw error if compensation does not fall within bounds
        if abs(f - 1) >= MAX_FCORR:
            raise PeakFindingFailed(
                "Compensation OOB",
                significance=sig,
                resolution=resolution,
                dt1=dt1,
                dt2=_dt1,
                dt=dt,
                df=f - 1,
            )

        # Stop if resolution met, otherwise refine resolution
        if resolution == target_resolution:
            break

        # Terminate immediately if 'quick' results desired
        if quick:
            break

        # Stop attempting frequency compensation if low enough
        if abs(df1) < TARGET_DF:
            do_frequency_compensation = False
            logger.debug("      Disabling frequency compensation.")

        # Update for next iteration
        bts = (bts - dt1) / (1 + df1)
        resolution /= separation_duration / duration / RES_REFINE_FACTOR
        resolution = max(resolution, target_resolution)
        curr_iteration += 1

    df = f - 1
    logger.debug("      Returning results.")
    return dt, df


@profile
def fpfind(
    alice: list,
    bob: list,
    num_wraps: int,
    num_bins: int,
    resolution: float,
    target_resolution: float,
    threshold: float,
    separation_duration: float,
    precompensations: list,
    precompensation_fullscan: bool = False,
    quick: bool = False,
    do_frequency_compensation: bool = True,
):
    """Performs fpfind procedure.

    This is effectively a wrapper to 'time_freq' that performs the actual
    frequency compensation routine, but also includes a frequency
    precompensation step, as well as potentially other further refinements
    (currently unimplemented).

    Timestamps must already by normalized to starting time of 0ns. Whether
    this starting time is also common reference between 'ats' and 'bts' is up to
    implementation.

    Args:
        alice: Reference timestamps, in 'a1' format.
        bob: Target timestamps, in 'a1' format.
        num_wraps: Number of cross-correlations to overlay, usually 1.
        num_bins: Numbers of bins to use in the FFT.
        resolution: Initial resolution of cross-correlation.
        target_resolution: Target resolution desired from routine.
        threshold: Height of peak to discriminate as signal, in units of dev.
        separation_duration: Separation of cross-correlations, in units of ns.
        precompensations: List of precompensations to apply.
        precompensation_fullscan: Specify whether to continue even after peak found.
    """
    df0s = precompensations

    # Go through all precompensations
    for df0 in df0s:
        logger.debug("  Applied initial %.4f ppm precompensation.", df0 * 1e6)

        # Apply frequency precompensation df0
        dt = 0
        f = 1 + df0
        try:
            dt1, df1 = time_freq(
                alice,
                (bob - dt) / f,
                num_wraps,
                num_bins,
                resolution,
                target_resolution,
                threshold,
                separation_duration,
                quick=quick,
                do_frequency_compensation=do_frequency_compensation,
            )
        except ValueError as e:
            logger.info(f"Peak finding failed, {df0*1e6:7.3f} ppm: {str(e)}")
            continue

        # Refine estimates, using the same recursive relations
        dt += f * dt1
        f *= 1 + df1
        logger.debug("  Applied another %.4f ppm compensation.", df1 * 1e6)
        e = PeakFindingFailed(
            "",
            resolution=target_resolution,
            dt=dt,
            df=f - 1,
        )
        logger.info(f"Peak found, precomp: {df0*1e6:7.3f} ppm: {str(e)}")
        # TODO: Justify the good enough frequency value
        # TODO(2024-01-31): Add looping code to customize refinement steps.
        if precompensation_fullscan:
            continue
        break

    # No appropriate frequency compensation found
    else:
        raise ValueError("No peak found!")  # TODO

    df = f - 1
    return dt, df


def generate_precompensations(start, stop, step, ordered=False) -> list:
    """Returns set of precompensations to apply before fpfind.

    The precompensations are in alternating positive/negative to allow
    scanning, e.g. for 10ppm, the sequence of values are:
    0ppm, 10ppm, -10ppm, 20ppm, -20ppm, ...

    Examples:
        >>> generate_precompensations(0, 0, 0)
        [0]
        >>> generate_precompensations(1, 10, 5)
        [1, 6, -4, 11, -9]
        >>> generate_precompensations(1, 1, 5)
        [1, 6, -4]
        >>> generate_precompensations(1, 10, 5, ordered=True)
        [1, 6, 11]
    """
    # Zero precompensation if invalid step supplied
    if step == 0:
        return [0]

    # Prepare pre-compensations
    if ordered:
        df0s = np.arange(0, int((stop - start) // step + 1))  # conventional
    else:
        df0s = np.arange(1, int(stop // step + 1) * 2) // 2  # flip-flop
        df0s[::2] *= -1

    df0s = df0s.astype(np.float64) * step
    df0s = df0s + start
    return df0s


# fmt: on
def main():
    global \
        _ENABLE_BREAKPOINT, \
        _DISABLE_DOUBLING, \
        _NUM_WRAPS_LIMIT, \
        TARGET_DF, \
        RES_REFINE_FACTOR
    script_name = Path(sys.argv[0]).name

    # Disable Black formatting
    # fmt: off

    def make_parser(help_verbosity=1):
        # Adapted from <https://discuss.python.org/t/advanced-help-for-argparse/20319/2>
        def adv(description):
            return description if help_verbosity >= 2 else configargparse.SUPPRESS
        def advv(description):
            return description if help_verbosity >= 3 else configargparse.SUPPRESS
        def advvv(description):
            return description if help_verbosity >= 4 else configargparse.SUPPRESS

        parser = configargparse.ArgumentParser(
            add_config_file_help=help_verbosity >= 2,
            default_config_files=[f"{script_name}.default.conf"],
            description=parse_docstring_description(__doc__),
            formatter_class=ArgparseCustomFormatter,
            add_help=False,
        )

        # Display arguments (group with defaults)
        pgroup_config = parser.add_argument_group("display/configuration")
        pgroup_config.add_argument(
            "-h", "--help", action="count", default=0,
            help="Show this help message, with incremental verbosity, e.g. up to -hhh")
        pgroup_config.add_argument(
            "-v", "--verbosity", action="count", default=0,
            help="Specify debug verbosity, e.g. -vv for more verbosity")
        pgroup_config.add_argument(
            "-L", "--logging", metavar="",
            help=adv("Log to file, if specified. Log level follows verbosity"))
        pgroup_config.add_argument(
            "--config", metavar="", is_config_file_arg=True,
            help=adv("Path to configuration file"))
        pgroup_config.add_argument(
            "--save", metavar="", is_write_out_config_file_arg=True,
            help=adv("Path to configuration file for saving, then immediately exit"))
        pgroup_config.add_argument(
            "--experiment", action="store_true",
            help=advvv("Enable debugging mode (needs 'python3 -im fpfind.fpfind')"))

        # Timestamp importing arguments
        pgroup_ts = parser.add_argument_group("importing timestamps")
        pgroup_ts.add_argument(
            "-t", "--reference", metavar="",
            help="Timestamp file in 'a1' format, from low-count side (reference)")
        pgroup_ts.add_argument(
            "-T", "--target", metavar="",
            help="Timestamp file in 'a1' format, from high-count side")
        pgroup_ts.add_argument(
            "-X", "--legacy", action="store_true",
            help="Parse raw timestamps in legacy mode (default: %(default)s)")
        pgroup_ts.add_argument(
            "-Z", "--skip-duration", metavar="", type=float, default=0,
            help=adv("Specify initial duration to skip, in seconds (default: %(default)s)"))

        # Epoch importing arguments
        pgroup_ep = parser.add_argument_group("importing epochs")
        pgroup_ep.add_argument(
            "-d", "--sendfiles", metavar="",
            help="SENDFILES, from low-count side (reference)")
        pgroup_ep.add_argument(
            "-D", "--t1files", metavar="",
            help="T1FILES, from high-count side")
        pgroup_ep.add_argument(
            "-e", "--first-epoch", metavar="",
            help=adv("Specify filename of first overlapping epoch, optional"))
        pgroup_ep.add_argument(
            "-z", "--skip-epochs", metavar="", type=int, default=0,
            help=adv("Specify number of initial epochs to skip (default: %(default)d)"))

        # Epoch importing arguments
        pgroup_chselect = parser.add_argument_group("channel selection")
        pgroup_chselect.add_argument(
            "-m", "--reference-pattern", metavar="", type=int,
            help=adv("Pattern mask for selecting detector events from low-count side"))
        pgroup_chselect.add_argument(
            "-M", "--target-pattern", metavar="", type=int,
            help=adv("Pattern mask for selecting detector events from high-count side"))

        # fpfind parameters
        pgroup_fpfind = parser.add_argument_group("fpfind parameters")
        pgroup_fpfind.add_argument(
            "--disable-comp", action="store_true",
            help=advv("Disables frequency compensation entirely"))
        pgroup_fpfind.add_argument(
            "-k", "--num-wraps", metavar="", type=int, default=1,
            help=adv("Specify number of arrays to wrap (default: %(default)d)"))
        pgroup_fpfind.add_argument(
            "--num-wraps-limit", metavar="", type=int, default=4,
            help=advv("Enables and specifies peak recovery 'num_wraps' limit (default: %(default)d)"))
        pgroup_fpfind.add_argument(
            "-q", "--buffer-order", metavar="", type=int, default=26,
            help="Specify FFT buffer order, N = 2**q (default: %(default)d)")
        pgroup_fpfind.add_argument(
            "-R", "--initial-res", metavar="", type=int, default=16,
            help="Specify initial coarse timing resolution, in units of ns (default: %(default)dns)")
        pgroup_fpfind.add_argument(
            "-r", "--final-res", metavar="", type=int, default=1,
            help=adv("Specify desired fine timing resolution, in units of ns (default: %(default)dns)"))
        pgroup_fpfind.add_argument(
            "-s", "--separation", metavar="", type=float, default=6,
            help=adv("Specify width of separation, in units of epochs (default: %(default).1f)"))
        pgroup_fpfind.add_argument(
            "-S", "--peak-threshold", metavar="", type=float, default=6,
            help=adv("Specify the statistical significance threshold (default: %(default).1f)"))
        pgroup_fpfind.add_argument(
            "--disable-doubling", action="store_true",
            help=advv("Disabling automatic resolution doubling during initial peak search"))
        pgroup_fpfind.add_argument(
            "--freq-threshold", metavar="", type=float, default=0.1,
            help=advv("Specify the threshold for frequency calculation, in units of ppb (default: %(default).1f)"))
        pgroup_fpfind.add_argument(
            "--convergence-rate", metavar="", type=float, default=np.sqrt(2),
            help=advv("Specify the rate of fpfind convergence (default: %(default).4f)"))
        pgroup_fpfind.add_argument(
            "-f", "--quick", action="store_true",
            help=advv("Returns the first iteration results immediately"))
        pgroup_fpfind.add_argument(
            "-V", "--output", metavar="", type=int, default=0, choices=range(1<<5),
            help=adv(f"{ArgparseCustomFormatter.RAW_INDICATOR}"
                "Specify output verbosity. Results are tab-delimited (default: %(default)d)\n"
                "- Setting bit 0 inverts the freq and time compensations\n"
                "- Setting bit 1 changes freq units, from abs to 2^-34\n"
                "- Setting bit 2 removes time compensation\n"
                "- Setting bit 3 changes time units, from 1ns to 1/8ns\n"
                "- Setting bit 4 adds first epoch used"
            )
        )

        # Timing pre-compensation parameters
        #
        # This is used for further fine-tuning of inter-bin timing values, or an overall shift
        # to perform a rough timing correction. Since this is a rarely used option, and to
        # avoid impacting runtime of the original algorithm, this timing precompensation is
        # performed *before* the shifts in frequency.
        #
        # To align with the compensation terminology elucidated in the output inversion
        # explanation, the reference is compensated for by the timing pre-compensation, i.e. alice.
        pgroup_tprecomp = parser.add_argument_group("timing precompensation")
        pgroup_tprecomp.add_argument(
            "--dt", metavar="", type=float, default=0,
            help=advvv("Specify initial timing shift, in units of ns (default: %(default)dns)"))
        pgroup_tprecomp.add_argument(
            "--dt-use-bins", action="store_true",
            help=advvv("Change dt units, from ns to timing resolution (i.e. -R)"))

        # Frequency pre-compensation parameters
        pgroup_precomp = parser.add_argument_group("frequency precompensation")
        pgroup_precomp.add_argument(
            "-P", "--precomp-enable", action="store_true",
            help="Enable precompensation scanning")
        pgroup_precomp.add_argument(
            "--df", "--precomp-start", metavar="", type=float, default=0.0,
            help=adv("Specify the precompensation value (default: 0ppm)"))
        pgroup_precomp.add_argument(
            "--precomp-step", metavar="", type=float, default=0.1e-6,
            help=adv("Specify the step value (default: 0.1ppm)"))
        pgroup_precomp.add_argument(
            "--precomp-stop", metavar="", type=float, default=20e-6,
            help=adv("Specify the max scan range, one-sided (default: 20ppm)"))
        pgroup_precomp.add_argument(
            "--precomp-ordered", action="store_true",
            help=advv("Test precompensations in increasing order (default: %(default)s)"))
        pgroup_precomp.add_argument(
            "--precomp-fullscan", action="store_true",
            help=advv("Force all precompensations to be tested (default: %(default)s)"))

        return parser

    # fmt: on
    # Parse arguments
    parser = make_parser()
    args = parser.parse_args()

    # Print advanced help if specified
    if args.help > 0:
        parser = make_parser(help_verbosity=args.help)
        parser.print_help(sys.stderr)
        sys.exit(1)

    # Check whether options have been supplied, and print help otherwise
    args_sources = parser.get_source_to_settings_dict().keys()
    config_supplied = any(map(lambda x: x.startswith("config_file"), args_sources))
    if len(sys.argv) == 1 and not config_supplied:
        parser.print_help(sys.stderr)
        sys.exit(1)

    # Set logging level and log arguments
    if args.logging is not None:
        set_logfile(logger, args.logging, human_readable=True)
    logger.setLevel(verbosity2level(args.verbosity))
    logger.info("%s", args)

    # Set experimental mode
    if args.experiment:
        _ENABLE_BREAKPOINT = True

    # Set timing recovery mode
    if args.num_wraps_limit > 0:
        _NUM_WRAPS_LIMIT = int(args.num_wraps_limit)

    # Set frequency threshold
    if args.freq_threshold > 0:
        TARGET_DF = args.freq_threshold * 1e-9

    # Set convergence rate
    if args.convergence_rate > 0:
        RES_REFINE_FACTOR = args.convergence_rate

    # Disable resolution doubling, if option specified
    if args.disable_doubling:
        _DISABLE_DOUBLING = True

    # Verify minimum duration has been imported
    num_bins = 1 << args.buffer_order
    Ta = args.initial_res * num_bins * args.num_wraps
    Ts = args.separation * Ta
    minimum_duration = (args.separation + 2) * Ta
    logger.debug(
        "Reading timestamps...",
        extra={
            "details": [
                f"Required duration: {minimum_duration*1e-9:.1f}s "
                f"(cross-corr {Ta*1e-9:.1f}s)",
            ]
        },
    )

    # fmt: off

    # Obtain timestamps needed for fpfind
    #   alice: low count side - chopper - HeadT2 - sendfiles (reference)
    #   bob: high count side - chopper2 - HeadT1 - t1files
    if args.sendfiles is not None and args.t1files is not None:
        logger.info("  Reading from epoch directories...")
        _is_reading_ts = False

        # +1 epoch specified for use as buffer for frequency compensation
        required_epochs = np.ceil(minimum_duration/EPOCH_LENGTH).astype(np.int32) + args.skip_epochs + 1

        # Automatically choose first overlapping epoch if not supplied manually
        first_epoch, available_epochs = get_first_overlapping_epoch(
            args.sendfiles, args.t1files,
            first_epoch=args.first_epoch, return_length=True,
        )
        logger.debug("  ", extra={"details": [
            f"Available: {available_epochs:d} epochs "
            f"(need {required_epochs:d})",
            f"First epoch: {first_epoch}",
        ]})
        if available_epochs < required_epochs:
            logger.warning("  Insufficient epochs")

        # Read epochs
        alice, aps = get_timestamp_pattern(
            args.sendfiles, "T2",
            first_epoch, args.skip_epochs, required_epochs - args.skip_epochs)
        bob, bps = get_timestamp_pattern(
            args.t1files, "T1",
            first_epoch, args.skip_epochs, required_epochs - args.skip_epochs)

    elif args.target is not None and args.reference is not None:
        logger.info("  Reading from timestamp files...")
        _is_reading_ts = True

        # Get only required events
        tas, tae = read_a1_start_end(args.reference, legacy=args.legacy)
        tbs, tbe = read_a1_start_end(args.target, legacy=args.legacy)
        start = max(tas, tbs)
        end = min(tae, tbe)
        required_duration = (minimum_duration + Ta) * 1e-9 + args.skip_duration  # seconds
        if end - start < required_duration :
            logger.warning("  Insufficient timestamp events")

        alice, aps = read_a1(args.reference, legacy=args.legacy, end=start+required_duration)
        bob, bps = read_a1(args.target, legacy=args.legacy, end=start+required_duration)

    else:
        logger.error("Timestamp files/epochs must be supplied with -tT/-dD")
        sys.exit(1)

    # Select events only from specified channels
    if args.reference_pattern is not None:
        alice = alice[(aps & args.reference_pattern).astype(bool)]
    if args.target_pattern is not None:
        bob = bob[(bps & args.target_pattern).astype(bool)]

    # Normalize timestamps to common time reference near start, so that
    # frequency compensation will not shift the timing difference too far
    skip = args.skip_duration if _is_reading_ts else 0
    alice, bob = normalize_timestamps(alice, bob, skip=skip)
    logger.debug(
        "  Read %d and %d events from reference and compensating side.",
        len(alice), len(bob), extra={"details": [
            "Reference timing range: "
            f"[{alice[0]*1e-9:.2f}, {alice[-1]*1e-9:.2f}]s",
            "Compensating timing range: "
            f"[{bob[0]*1e-9:.2f}, {bob[-1]*1e-9:.2f}]s",
            f"(skipped {skip:.2f}s)",
        ]
    })

    # Prepare frequency pre-compensations
    precompensations = [args.df]
    if args.precomp_enable:
        logger.debug("Generating frequency precompensations...")
        precompensations = generate_precompensations(
            args.df,
            args.precomp_stop,
            args.precomp_step,
            ordered=args.precomp_ordered,
        )
        logger.debug(
            "  Prepared %d precompensation(s): %s... ppm",
            len(precompensations),
            ",".join(map(lambda p: f"{p*1e6:g}", precompensations[:3])),
        )

    # Perform timing pre-compensation (experimental command)
    # See notes in CLI definition.
    if args.dt != 0:
        alice = alice + args.dt * (args.initial_res if args.dt_use_bins else 1)

    # Start fpfind
    logger.debug("Running fpfind...")
    dt, df = fpfind(
        alice, bob,
        num_wraps=args.num_wraps,
        num_bins=num_bins,
        resolution=args.initial_res,
        target_resolution=args.final_res,
        threshold=args.peak_threshold,
        separation_duration=Ts,
        precompensations=precompensations,
        precompensation_fullscan=args.precomp_fullscan,
        quick=args.quick,
        do_frequency_compensation=not args.disable_comp,
    )

    # To understand the options below, we first clarify some definitions:
    #   - To apply a frequency difference is to apply an additional clock
    #     skew to the target timestamp events, i.e.  ts * (1 + df).
    #       - In the context of timing difference:   ts + dt
    #   - To compensate for a frequency difference is to undo the clock
    #     skew present in the target events, i.e.    ts / (1 + df)
    #       - In the context of timing difference:   ts - dt
    #   - 'fpfind' calculates the *compensation* frequency for the target,
    #     and the *compensation* timing prior to frequency compensation,
    #     i.e.  (ts - dt) / (1 + df)
    #       - This is equivalent to applying the frequency difference to
    #         the reference, *then* applying the timing difference,
    #         i.e.  Ts * (1 + df) + dt
    #   - To invert this relation is to have the results be *applied* to
    #     the target (vis a vis *compensated* for the reference), i.e.
    #           ts * (1 + dF) + dT
    #           (Ts - dT) / (1 + dF)
    #       - This convention is followed by the qcrypto stack, where
    #         'freqcd' first applies the freq diff to the target, before
    #         downstream 'costream' corrects and tracks the timing diff.
    #   - This implies the inversion is represented by:
    #        dT = -dt / (1 + df)
    #        dF = 1 / (1 + df) - 1

    # Modify output based on flags
    flag = args.output
    if flag & 0b0001:
        dt = -dt / (1 + df)
        df = 1 / (1 + df) - 1

    # Convert into special bases
    if flag & 0b0010:
        df = f"{round(df * (1 << 34)):d}"  # units of 2**-34 as per freqcd
    if flag & 0b1000:
        dt *= 8  # units of 1/8ns as per qcrypto

    # Generate output, and add timing value if requested
    output = f"{df}\t"
    if not (flag & 0b0100):
        output += f"{round(dt):d}\t"
    if (flag & 0b10000) and not _is_reading_ts:
        output += f"{first_epoch}\t"
    output = output.rstrip()
    print(output, file=sys.stdout)  # newline auto-added


if __name__ == "__main__":
    main()
