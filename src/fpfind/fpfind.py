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

import sys
from pathlib import Path

import configargparse
import numpy as np
import numpy.typing as npt
from typing_extensions import TypeAlias

import fpfind.lib._logging as logging
from fpfind import VERSION
from fpfind.lib.constants import (
    EPOCH_LENGTH,
    MAX_FCORR,
    MAX_TIMING_RESOLUTION_NS,
    NTP_MAXDELAY_NS,
    FrequencyCompensation,
    PeakFindingFailed,
)
from fpfind.lib.parse_timestamps import read_a1, read_a1_start_end
from fpfind.lib.utils import (
    ArgparseCustomFormatter,
    fold_histogram,
    get_first_overlapping_epoch,
    get_statistics,
    get_timestamp_pattern,
    get_timing_delay_fft,
    histogram_fft3,
    match_dts,
    normalize_timestamps,
    parse_docstring_description,
    round,
)

logger, log = logging.get_logger("fpfind")

# Disables resolution doubling during initial peak search, used to converge
# on specific fpfind parameters for peak search. This procedure is a relatively
# cheap measure to increase peak finding yield, and thus is unlikely needed
# in production.
# For internal use only.
_DISABLE_DOUBLING = False

# Toggles interruptible FFT
ENABLE_INTERRUPT = False

# Type aliases
NDArrayFloat: TypeAlias = npt.NDArray[np.floating]
NDArrayNumber: TypeAlias = npt.NDArray[np.number]


# Main algorithm
def time_freq(
    ats: NDArrayNumber,
    bts: NDArrayNumber,
    k0: int,
    N: int,
    r0: float,
    r_target: float,
    S0: float,
    max_dt: float,
    max_df: float,
    df_target: float,
    Ts: float,
    convergence_rate: float,
    quick: bool,
    do_frequency_compensation: FrequencyCompensation = FrequencyCompensation.ENABLE,
):
    """Perform the actual frequency compensation routine.

    Timestamps must already by normalized to starting time of 0ns. Whether
    this starting time is also common reference between 'ats' and 'bts' is up to
    implementation.

    Note some differences with original implementation: the 'separation_duration'
    is left constant.

    Args:
        ats: Timestamps of reference side, in units of ns.
        bts: Timestamps of compensating side, in units of ns.
        k0: Number of cross-correlations to overlay, usually 1.
        N: Numbers of bins to use in the FFT.
        r0: Initial resolution of cross-correlation.
        r_target: Target resolution desired from routine.
        S0: Height of peak to discriminate as signal, in units of dev.
        Ts: Separation of cross-correlations, in units of ns.
    """
    # Optional new behaviour where the significance evaluation is deferred,
    # instead looking for coincidences in timing differences
    perform_liberal_match = S0 == 0
    r = _r0 = r0
    k = k0
    Ta = r0 * N * k0

    log(1).debug("Performing peak searching...")
    log(2).debug(
        "Parameters:",
        f"Bins: 2^{np.int32(np.log2(N)):d}",
        f"Bin width: {r0:.0f}ns",
        f"Number of wraps: {k0:d}",
        f"Target resolution: {r_target}ns",
    )

    # Refinement loop, note resolution/duration will change during loop
    dt = 0
    f = 1
    iter = 1
    perform_coarse_finding = True
    prev_dt1 = 0  # cached dt1 (not _dt1!)

    while True:
        # Use previously cached base resolution
        if perform_liberal_match:
            r = _r0

        # Dynamically adjust 'k' to avoid event over- and under-flow
        k_max = int(np.floor(Ta / (r * N)))
        k_act = min(k, max(k_max, 1))  # required because 'k_max' may < 1
        Ta_act = k_act * r * N
        Ta_act = min(Ta_act, Ta)
        log(2).debug(f"Iteration {iter} (r={r:.1f}ns, k={k_act:d})")

        # Perform cross-correlation
        log(3).debug(f"Performing earlier xcorr (range: [0.00, {Ta_act * 1e-9:.2f}]s)")
        xs, ys = histogram_fft3(
            ats, bts, 0, Ta_act, N, r, Ta, interruptible=ENABLE_INTERRUPT
        )

        # Calculate timing delay
        dt1 = get_timing_delay_fft(ys, xs)[0]  # get smaller candidate
        sig = get_statistics(ys, r).significance

        # Dynamic search for resolution
        dt1s_early = []
        while perform_liberal_match or perform_coarse_finding:
            log(4).debug(
                f"Peak: S = {sig:.3f}, dt = {dt1}ns (resolution = {r:.0f}ns)",
            )

            # Deviation zero, due flat cross-correlation
            # Hint: Increase number of bins
            if sig == 0:
                raise PeakFindingFailed(
                    "Bin saturation  ",
                    significance=sig,
                    resolution=r,
                    dt1=dt1,
                )

            if abs(dt1) < max_dt:
                dt1s_early.append((dt1, r, sig))

            # Check if peak threshold exceeded
            if not perform_liberal_match:
                if sig >= S0:
                    log(5).debug("Accepted")
                    break
                log(5).debug("Rejected")

            if _DISABLE_DOUBLING:
                raise PeakFindingFailed(
                    "Low significance",
                    significance=sig,
                    resolution=r,
                    dt1=dt1,
                )

            # If peak rejected, merge contiguous bins to double
            # the resolution of the peak search
            r *= 2
            xs, ys = fold_histogram(xs, ys, 2)
            dt1 = get_timing_delay_fft(ys, xs)[0]
            sig = get_statistics(ys, r).significance

            # Catch runaway resolution doubling, limited by
            if r > MAX_TIMING_RESOLUTION_NS:
                if perform_liberal_match:
                    if len(dt1s_early) == 0:
                        raise PeakFindingFailed(
                            "Time delay OOB  ",
                            significance=sig,
                            resolution=r,
                            dt1=dt1,
                        )
                    break
                raise PeakFindingFailed(
                    "Resolution OOB  ",
                    significance=sig,
                    resolution=r,
                    dt1=dt1,
                )

        # Calculate some thresholds to catch when peak was likely not found
        buffer = 1
        threshold_dt = buffer * max(abs(prev_dt1), 1)

        # If the duration has too few coincidences, the peak may not
        # show up at all. A peak is likely to have already been found, if
        # the current iteration is more than 1. Attempt to retry with
        # larger 'num_wraps' if enabled in settings.
        if (
            not perform_liberal_match
            and not perform_coarse_finding
            and abs(dt1) > threshold_dt
        ):
            log(3).warning(
                "Interrupted due spurious signal:",
                f"early dt       = {dt1:10.0f} ns",
                f"threshold dt   = {threshold_dt:10.0f} ns",
            )
            k *= 2
            if k <= k_max:
                log(3).debug(
                    f"Reattempting with k = {k:d} due missing peak.",
                )
                continue
            break

        # Catch if timing delay exceeded
        if not perform_liberal_match and abs(dt1) > max_dt:
            raise PeakFindingFailed(
                "Time delay OOB  ",
                significance=sig,
                resolution=r,
                dt1=dt1,
            )

        _dt1 = dt1
        df1 = 0  # default values
        if perform_liberal_match or (
            do_frequency_compensation is not FrequencyCompensation.DISABLE
        ):
            # Use previously cached base resolution
            if perform_liberal_match:
                r = _r0

            # Perform cross-correlation
            log(3).debug(
                f"Performing later xcorr (range: [{Ts * 1e-9:.2f}, {(Ts + Ta_act) * 1e-9:.2f}]s)",
            )
            xs, _ys = histogram_fft3(
                ats, bts, Ts, Ta_act, N, r, Ta, interruptible=ENABLE_INTERRUPT
            )

            # Calculate timing delay for late set of timestamps
            _dt1 = get_timing_delay_fft(_ys, xs)[0]
            sig = get_statistics(_ys, r).significance

            # Attempt similar search for late timing difference
            dt1s_late = []
            while perform_liberal_match:
                log(4).debug(
                    f"Peak: S = {sig:.3f}, dt = {_dt1}ns (resolution = {r:.0f}ns)",
                )
                if abs(_dt1) < max_dt:
                    dt1s_late.append((_dt1, r, sig))

                r *= 2
                xs, _ys = fold_histogram(xs, _ys, 2)
                _dt1 = get_timing_delay_fft(_ys, xs)[0]
                sig = get_statistics(_ys, r).significance

                # Check intersections here
                if r > MAX_TIMING_RESOLUTION_NS:
                    if len(dt1s_late) == 0:
                        raise PeakFindingFailed(
                            "Time delay OOB  ",
                            significance=sig,
                            resolution=r,
                            dt1=_dt1,
                        )

                    allowed_dt_diff = Ts * max_df * 1e9
                    dt1, _dt1, r = match_dts(dt1s_early, dt1s_late, allowed_dt_diff)
                    break

            # Evaluate measured frequency difference
            df1 = (_dt1 - dt1) / Ts

            # Some guard rails to make sure results make sense
            # Something went wrong with peak searching, to return intermediate
            # results which are likely near correct values.
            threshold_df = buffer * max(abs(f - 1), 1e-9)
            if f == 1:
                threshold_df = buffer * MAX_FCORR * 1e6  # [ppm]

            if (
                not perform_liberal_match
                and not perform_coarse_finding
                and (abs(_dt1) > threshold_dt or abs(df1) > threshold_df)
            ):
                log(3).warning(
                    "Interrupted due spurious signal:",
                    f"early dt       = {dt1:10.0f} ns",
                    f"late dt        = {_dt1:10.0f} ns",
                    f"threshold dt   = {threshold_dt:10.0f} ns",
                    f"current df     = {df1 * 1e6:10.4f} ppm",
                    f"threshold df   = {threshold_df * 1e6:10.4f} ppm",
                )
                k *= 2
                if k <= k_max:
                    log(3).debug(
                        f"Reattempting with k = {k:d} due missing peak.",
                    )
                    continue

                break  # terminate if recovery not enabled

            # Catch if timing delay exceeded
            if not perform_liberal_match and abs(_dt1) > max_dt:
                raise PeakFindingFailed(
                    "Time delay OOB  ",
                    significance=sig,
                    resolution=r,
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
        log(3).debug(
            "Calculated timing delays:",
            f"early dt       = {dt1:10.0f} ns",
            f"late dt        = {_dt1:10.0f} ns",
            f"accumulated dt = {dt:10.0f} ns",
            f"current df     = {df1 * 1e6:10.4f} ppm",
            f"accumulated df = {(f - 1) * 1e6:10.4f} ppm",
        )

        # Throw error if compensation does not fall within bounds
        if abs(f - 1) >= max_df:
            raise PeakFindingFailed(
                "Compensation OOB",
                significance=sig,
                resolution=r,
                dt1=dt1,
                dt2=_dt1,
                dt=dt,
                df=f - 1,
            )

        # Stop if resolution met, otherwise refine resolution
        if r <= r_target:
            break

        # Stop liberal search if current resolution is expected to
        # return a reasonable result
        if r == _r0:
            perform_liberal_match = False

        # Terminate immediately if 'quick' results desired
        if quick:
            break

        # Stop attempting frequency compensation if low enough
        if abs(df1) < df_target and (
            do_frequency_compensation is FrequencyCompensation.ENABLE
        ):
            do_frequency_compensation = FrequencyCompensation.DISABLE
            log(3).debug("Disabling frequency compensation.")

        # if perform_coarse_finding:
        #     N //= int(np.round(r / r0))

        # Update for next iteration
        max_dt1 = max(abs(dt1), abs(_dt1))
        if max_dt1 != 0:
            prev_dt1 = max_dt1
        bts = (bts - dt1) / (1 + df1)
        r = max(r / convergence_rate, r_target)
        iter += 1
        perform_coarse_finding = False  # disable coarse search, i.e. iter > 1

        # Cache resolution for the next iteration
        if perform_liberal_match:
            _r0 = r

    df = f - 1
    log(3).debug("Returning results.")
    return dt, df


def fpfind(
    alice: NDArrayNumber,
    bob: NDArrayNumber,
    num_wraps: int,
    num_bins: int,
    resolution: float,
    target_resolution: float,
    threshold: float,
    max_dt: float,
    max_df: float,
    df_target: float,
    separation_duration: float,
    convergence_rate: float,
    precompensations: list,
    precompensation_fullscan: bool = False,
    quick: bool = False,
    do_frequency_compensation: FrequencyCompensation = FrequencyCompensation.ENABLE,
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
        log(1).debug(
            f"Applied initial {df0 * 1e6:.4f} ppm precompensation.",
        )

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
                max_dt,
                max_df,
                df_target,
                separation_duration,
                quick=quick,
                convergence_rate=convergence_rate,
                do_frequency_compensation=do_frequency_compensation,
            )
        except ValueError as e:
            log(0).info(f"Peak finding failed, {df0 * 1e6:7.3f} ppm: {str(e)}")
            continue

        # Refine estimates, using the same recursive relations
        dt += f * dt1
        f *= 1 + df1
        log(1).debug(
            f"Applied another {df1 * 1e6:.4f} ppm compensation.",
        )
        e = PeakFindingFailed(
            "",
            resolution=target_resolution,
            dt=dt,
            df=f - 1,
        )
        log(0).info(f"Peak found, precomp: {df0 * 1e6:7.3f} ppm: {str(e)}")
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
    global ENABLE_INTERRUPT, _DISABLE_DOUBLING
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
        pgroup = parser.add_argument_group("display/configuration")
        pgroup.add_argument(
            "-h", "--help", action="count", default=0,
            help="Show this help message, with incremental verbosity, e.g. up to -hhh")
        pgroup.add_argument(
            "-v", "--verbosity", action="count", default=0,  # retained for backward compatibility
            help=configargparse.SUPPRESS)  # black hole for deprecated option
        pgroup.add_argument(
            "-p", "--quiet", action="count", default=0,
            help="Reduce log verbosity, e.g. -pp for less verbosity")
        pgroup.add_argument(
            "-L", "--logging", metavar="",
            help=adv("Log to file, if specified. Log level follows verbosity"))
        pgroup.add_argument(
            "--config", metavar="", is_config_file_arg=True,
            help=adv("Path to configuration file"))
        pgroup.add_argument(
            "--save", metavar="", is_write_out_config_file_arg=True,
            help=adv("Path to configuration file for saving, then immediately exit"))
        pgroup.add_argument(
            "-I", "--interruptible", action="store_true",
            help=advv("Allow fpfind routine to be interrupted via SIGINT"))
        pgroup.add_argument(
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

        # Timestamp importing arguments
        pgroup = parser.add_argument_group("importing timestamps")
        pgroup.add_argument(
            "-t", "--reference", metavar="",
            help="Timestamp file in 'a1' format, from low-count side (reference)")
        pgroup.add_argument(
            "-T", "--target", metavar="",
            help="Timestamp file in 'a1' format, from high-count side")
        pgroup.add_argument(
            "-X", "--legacy", action="store_true",
            help="Parse raw timestamps in legacy mode (default: %(default)s)")
        pgroup.add_argument(
            "-Z", "--skip-duration", metavar="", type=float, default=0,
            help=adv("Specify initial duration to skip, in seconds (default: %(default)s)"))

        # Epoch importing arguments
        pgroup = parser.add_argument_group("importing epochs")
        pgroup.add_argument(
            "-d", "--sendfiles", metavar="",
            help="SENDFILES, from low-count side (reference)")
        pgroup.add_argument(
            "-D", "--t1files", metavar="",
            help="T1FILES, from high-count side")
        pgroup.add_argument(
            "-e", "--first-epoch", metavar="",
            help=adv("Specify filename of first overlapping epoch, optional"))
        pgroup.add_argument(
            "-z", "--skip-epochs", metavar="", type=int, default=0,
            help=adv("Specify number of initial epochs to skip (default: %(default)d)"))

        # Channel selection
        pgroup = parser.add_argument_group("channel selection")
        pgroup.add_argument(
            "-m", "--reference-pattern", metavar="", type=int,
            help=adv("Pattern mask for selecting detector events from low-count side"))
        pgroup.add_argument(
            "-M", "--target-pattern", metavar="", type=int,
            help=adv("Pattern mask for selecting detector events from high-count side"))

        # Timing compensation (pfind) parameters
        pgroup = parser.add_argument_group("timing compensation")
        pgroup.add_argument(
            "-k", "--num-wraps", metavar="", type=int, default=1,
            help=adv("Specify number of arrays to wrap (default: %(default)d)"))
        pgroup.add_argument(
            "-q", "--buffer-order", metavar="", type=int, default=26,
            help="Specify FFT buffer order, N = 2**q (default: %(default)d)")
        pgroup.add_argument(
            "-R", "--initial-res", metavar="", type=int, default=16,
            help="Specify initial coarse timing resolution, in units of ns (default: %(default)dns)")
        pgroup.add_argument(
            "-r", "--final-res", metavar="", type=int, default=1,
            help=adv("Specify desired fine timing resolution, in units of ns (default: %(default)dns)"))
        pgroup.add_argument(
            "--max-dt", metavar="", type=float, default=NTP_MAXDELAY_NS,
            help=advv("Expected maximum timing difference, in units of ns (default: 200ms)"))
        pgroup.add_argument(
            "-S", "--peak-threshold", metavar="", type=float, default=6,
            help=adv("Specify the statistical significance threshold (default: %(default).1f)"))
        pgroup.add_argument(
            "--disable-doubling", action="store_true",
            help=advv("Disabling automatic resolution doubling during initial peak search"))
        pgroup.add_argument(
            "--convergence-rate", metavar="", type=float,
            help=configargparse.SUPPRESS)  # black hole for deprecated option
        pgroup.add_argument(
            "-Q", "--convergence-order", metavar="", type=float, default=4,
            help=adv("Specify the reduction factor in timing uncertainty between iterations, larger = faster (default: %(default).4f)"))
        pgroup.add_argument(
            "-f", "--quick", action="store_true",
            help=advv("Returns the first iteration results immediately"))

        # Frequency compensation parameters
        pgroup = parser.add_argument_group("frequency compensation")
        pgroup.add_argument(
            "-s", "--separation", metavar="", type=float, default=6,
            help=adv("Specify width of separation, in units of epochs (default: %(default).1f)"))
        pgroup.add_argument(
            "--force-comp", action="store_true",
            help=advv("Forces frequency compensation even when no drift detected"))
        pgroup.add_argument(
            "--disable-comp", action="store_true",
            help=advv("Disables frequency compensation entirely"))
        pgroup.add_argument(
            "--max-df", metavar="", type=float, default=MAX_FCORR,
            help=advv("Expected maximum frequency difference (default: 122ppm)"))
        pgroup.add_argument(
            "--freq-threshold", metavar="", type=float, default=0.1,
            help=advv("Threshold for frequency calculation, in units of ppb (default: %(default).1f)"))

        # Timing pre-compensation parameters
        #
        # This is used for further fine-tuning of inter-bin timing values, or an overall shift
        # to perform a rough timing correction. Since this is a rarely used option, and to
        # avoid impacting runtime of the original algorithm, this timing precompensation is
        # performed *before* the shifts in frequency.
        #
        # To align with the compensation terminology elucidated in the output inversion
        # explanation, the reference is compensated for by the timing pre-compensation, i.e. alice.
        pgroup = parser.add_argument_group("timing precompensation")
        pgroup.add_argument(
            "--dt", metavar="", type=float, default=0,
            help=advv("Initial timing shift, in units of ns (default: %(default)dns)"))
        pgroup.add_argument(
            "--dt-use-bins", action="store_true",
            help=advv("Change dt units, from ns to timing resolution (i.e. -R)"))

        # Frequency pre-compensation parameters
        pgroup = parser.add_argument_group("frequency precompensation")
        pgroup.add_argument(
            "-P", "--precomp-enable", action="store_true",
            help="Enable precompensation scanning")
        pgroup.add_argument(
            "--df", "--precomp-start", metavar="", type=float, default=0.0,
            help=adv("Specify the precompensation value (default: 0ppm)"))
        pgroup.add_argument(
            "--precomp-step", metavar="", type=float, default=0.1e-6,
            help=adv("Specify the step value (default: 0.1ppm)"))
        pgroup.add_argument(
            "--precomp-stop", metavar="", type=float, default=20e-6,
            help=adv("Specify the max scan range, one-sided (default: 20ppm)"))
        pgroup.add_argument(
            "--precomp-ordered", action="store_true",
            help=advv("Test precompensations in increasing order (default: %(default)s)"))
        pgroup.add_argument(
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

    # Set log arguments
    if args.logging is not None:
        logging.set_logfile(logger, args.logging)

    # Set logging level
    # Default level should be DEBUG (instead of WARNING) for better
    # accessibility to the tool. '--verbosity' kept for legacy reasons.
    if args.quiet > 0:
        verbosity = max(2 - args.quiet, 0)
    elif args.verbosity == 0:
        verbosity = 2
    else:
        verbosity = args.verbosity
    logging.set_verbosity(logger, verbosity)
    log(0).info(f"fpfind {VERSION}")
    log(0).info(f"{args}")
    if args.quiet == 0 and args.verbosity > 0:
        log(0).warning(
            "'-v'/'--verbosity' has been deprecated, use "
            "'-p'/'--quiet' instead with inverted behaviour."
        )

    # Allow interrupting fpfind routine
    if args.interruptible:
        ENABLE_INTERRUPT = True

    # Set lower threshold limit to frequency compensation before disabling
    df_target = 1e-10
    if args.freq_threshold > 0:
        df_target = args.freq_threshold * 1e-9

    # Set convergence rate
    if args.convergence_rate is not None:
        log(0).error(
            "'--convergence-rate' has been removed, use "
            "'-Q'/'--convergence-order' instead with new behaviour."
        )
        sys.exit(1)

    # Enforce only one frequency compensation setting
    if args.disable_comp and args.force_comp:
        log(0).error("Only one of '--disable-comp' and '--force-comp' allowed.")
        sys.exit(1)
    do_frequency_compensation = FrequencyCompensation.ENABLE
    if args.disable_comp:
        do_frequency_compensation = FrequencyCompensation.DISABLE
    if args.force_comp:
        do_frequency_compensation = FrequencyCompensation.FORCE

    # Disable resolution doubling, if option specified
    if args.disable_doubling:
        _DISABLE_DOUBLING = True

    # Verify minimum duration has been imported
    num_bins = 1 << args.buffer_order
    Ta = args.initial_res * num_bins * args.num_wraps
    Ts = args.separation * Ta
    minimum_duration = (args.separation + 2) * Ta  # TODO: Check if this should be 1
    log(0).info(
        "Reading timestamps...",
        f"Required duration: {minimum_duration * 1e-9:.1f}s "
        f"(cross-corr {Ta * 1e-9:.1f}s)",
    )

    # fmt: off

    # Obtain timestamps needed for fpfind
    #   alice: low count side - chopper - HeadT2 - sendfiles (reference)
    #   bob: high count side - chopper2 - HeadT1 - t1files
    if args.sendfiles is not None and args.t1files is not None:
        log(1).info("Reading from epoch directories...")
        _is_reading_ts = False

        # +1 epoch specified for use as buffer for frequency compensation
        required_epochs = np.ceil(minimum_duration/EPOCH_LENGTH).astype(np.int32) + args.skip_epochs + 1

        # Automatically choose first overlapping epoch if not supplied manually
        first_epoch, available_epochs = get_first_overlapping_epoch(
            args.sendfiles, args.t1files,
            first_epoch=args.first_epoch, return_length=True,
        )  # type: ignore
        log(1).debug(
            "",
            f"Available: {available_epochs:d} epochs "
            f"(need {required_epochs:d})",
            f"First epoch: {first_epoch}",
        )
        if available_epochs < required_epochs:
            log(1).warning("Insufficient epochs")

        # Read epochs
        alice, aps = get_timestamp_pattern(
            args.sendfiles, "T2",
            first_epoch, args.skip_epochs, required_epochs - args.skip_epochs)
        bob, bps = get_timestamp_pattern(
            args.t1files, "T1",
            first_epoch, args.skip_epochs, required_epochs - args.skip_epochs)

    elif args.target is not None and args.reference is not None:
        log(1).info("Reading from timestamp files...")
        _is_reading_ts = True
        first_epoch = None

        # Get only required events
        tas, tae = read_a1_start_end(args.reference, legacy=args.legacy)
        tbs, tbe = read_a1_start_end(args.target, legacy=args.legacy)
        start = max(tas, tbs)
        end = min(tae, tbe)
        required_duration = (minimum_duration + Ta) * 1e-9 + args.skip_duration  # seconds
        if end - start < required_duration :
            log(1).warning("Insufficient timestamp events")

        alice, aps = read_a1(args.reference, legacy=args.legacy, end=start+required_duration)
        bob, bps = read_a1(args.target, legacy=args.legacy, end=start+required_duration)

    else:
        log(0).error("Timestamp files/epochs must be supplied with -tT/-dD")
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
    log(1).debug(
        f"Read {len(alice):d} and {len(bob):d} events from reference and compensating side.",
        "Reference timing range: "
        f"[{alice[0]*1e-9:.2f}, {alice[-1]*1e-9:.2f}]s",
        "Compensating timing range: "
        f"[{bob[0]*1e-9:.2f}, {bob[-1]*1e-9:.2f}]s",
        f"(skipped {skip:.2f}s)",
    )

    # Prepare frequency pre-compensations
    precompensations = [args.df]
    if args.precomp_enable:
        log(0).info("Generating frequency precompensations...")
        precompensations = generate_precompensations(
            args.df,
            args.precomp_stop,
            args.precomp_step,
            ordered=args.precomp_ordered,
        )
        text = ",".join(map(lambda p: f"{p*1e6:g}", precompensations[:3]))
        log(1).debug(
            f"Prepared {len(precompensations):d} precompensation(s): {text}... ppm",
        )

    # Perform timing pre-compensation (experimental command)
    # See notes in CLI definition.
    if args.dt != 0:
        alice = alice + args.dt * (args.initial_res if args.dt_use_bins else 1)

    # Start fpfind
    log(0).info("Running fpfind...")
    dt, df = fpfind(
        alice, bob,
        num_wraps=args.num_wraps,
        num_bins=num_bins,
        resolution=args.initial_res,
        target_resolution=args.final_res,
        threshold=args.peak_threshold,
        max_dt=args.max_dt,
        max_df=args.max_df,
        df_target=df_target,
        separation_duration=Ts,
        convergence_rate=args.convergence_order,
        precompensations=precompensations,
        precompensation_fullscan=args.precomp_fullscan,
        quick=args.quick,
        do_frequency_compensation=do_frequency_compensation,
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
