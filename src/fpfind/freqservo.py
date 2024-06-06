#!/usr/bin/env python3
"""Reads costream output and provides frequency compensation updates.

It reads the costream logging output (either written to '-n' argument, or to stdout otherwise),
and performs a linear estimate of the timing drift (frequency offset). These frequency offsets are compatible
with the formats read by 'freqcd'. See the documentation for 'freqservo:evaluate_freqshift' on how to
tune this estimation. Example scenario:

"freqcd -F freq.pipe -f freqdiff -o raws &
chopper2 -i raws -D t1dir -U &
costream -D t1dir -n log.pipe -V 5 &
freqservo -n log.pipe -V 5 -F freq.pipe -f freqdiff &".

"""

import logging
import sys
import time
from pathlib import Path

import configargparse
import numpy as np

from fpfind.lib import parse_epochs as eparser
from fpfind.lib.constants import EPOCH_DURATION
from fpfind.lib.utils import parse_docstring_description, ArgparseCustomFormatter

_LOGGING_FMT = "{asctime}\t{levelname:<7s}\t{funcName}:{lineno}\t| {message}"
logger = logging.getLogger(__name__)
if not logger.handlers:
    handler = logging.StreamHandler(stream=sys.stderr)
    handler.setFormatter(
        logging.Formatter(fmt=_LOGGING_FMT, datefmt="%Y%m%d_%H%M%S", style="{")
    )
    logger.addHandler(handler)
    logger.propagate = False

# Conversion factors for frequency units, i.e. 2^-34 <-> 1e-10
FCORR_DTOB = (1 << 34) / 10_000_000_000
FCORR_BTOD = 1/FCORR_DTOB

# Container to hold previously servoed time differences
DT_HISTORY = []

def parse(line, format=None):
    """Reads costream stdout and returns epoch and servoed time.

    The returned time is in units of ns. If 'format' is not supplied, 'line' is
    expected to be of the form '<EPOCH:str> <DT:float>\n', rather than a costream
    format.

    Default format (space-delimited):
        bead38d4 -77845

    costream V=4 format (space-delimited):
        epoch: bead378b, 2-evnts: 86180, 4-evnts: 2997, new bw4: 7, ft: 346, acc: 2876, true: 2997, 1-events: 86324

    costream V=5 format (tab-delimited):
        Epoch     Raw     Sent  Cmpr  Track   Acc   Coinc  Singles  (headers, not in output)
        bead38d4  116271  4144  7     -77845  3897  4144   111068
    """
    if format is None:
        epoch, dt = line.strip().split(" ")
        dt = float(dt)
    elif format == 5:
        tokens = line.split("\t")
        epoch = tokens[0]
        dt = int(tokens[4]) / 8  # convert 1/8ns -> 1ns
    elif format == 4:
        tokens = line.split(" ")
        epoch = tokens[1][:-1]
        dt = int(tokens[10][:-1]) / 8
    else:
        raise ValueError(f"Unrecognised costream format: '{format}'")
    logger.debug("Successfully parsed epoch = %s, dt = %.3f ns", epoch, dt)
    return epoch, dt

def evaluate_freqshift(epoch, dt, ignore=5, average=3, separation=12, cap=100e-9):
    """Returns desired clock skew correction, in absolute units.

    This function accumulates the dt values, and performs some degree of
    averaging (low-pass) to smoothen the frequency correction adjustment.

    Example:
        Suppose the dt history contains the set of time differences (in ns),
        with parameters (ignore=5, average=3, separation=12, cap=10e-9):

        dt_history = (
            { 43.250,  99.250,  148.625,  176.500, 197.375,}  # ignored
            {220.250, 229.750,  247.750,} 255.125, 268.000,   # initial sample
             271.875, 277.625,  281.625,  278.125, 282.250,
             292.625, 300.000, {304.000,  315.000, 314.125,}  # final sample
        )

        The averaged servoed timediff in the initial and final samples are
        232.58333 and 311.04166 ns, which in span of 12 epochs corresponds to
        a frequency drift of 12.178 ppb. The compensation frequency is thus
        -12.2 ppb, to the nearest ppb, and then capped to -10.0 ppb.

    Note:
        Assuming a bimodal clock skew, the algorithm used should not be a
        sliding mean, since the coincidence matching requires an accurate
        frequency compensation value. We instead collect a set of samples,
        then calculate the frequency difference.

        Race condition may be possible - no guarantees on the continuity
        of epochs. This is resolved with an epoch continuity check at every
        function call.

        This function was ported over from QKDServer.S15qkd.readevents[1].

    References:
        [1]: <https://github.com/s-fifteen-instruments/QKDServer/blob/master/S15qkd/readevents.py>
    """
    # Verify epochs are contiguous
    # Needed because costream may terminate prematurely, and no calls to
    # flush the dt history is made, resulting in large time gaps.
    if len(DT_HISTORY) > 0:
        # DT_HISTORY format: [(epoch:str, dt:float),...]
        prev_epoch = DT_HISTORY[-1][0]
        if eparser.epoch2int(prev_epoch) + 1 != eparser.epoch2int(epoch):
            DT_HISTORY.clear()

    # Collect epochs
    DT_HISTORY.append((epoch, dt))
    required = ignore + average + separation
    if len(DT_HISTORY) < required:
        return

    # Ignore first few epochs to allow freqcorr to propagate
    history = DT_HISTORY[ignore:required]  # make a copy
    epochs, dts = list(zip(*history))
    logger.debug("Active frequency correction triggered, measurements: %s", dts)

    # Calculate averaged frequency difference
    dt_early = np.mean(dts[:average])
    dt_late = np.mean(dts[-average:])
    dt_change = (dt_late - dt_early) * 1e-9  # convert to seconds
    df = dt_change / (separation * EPOCH_DURATION)
    df_toapply = 1/(1 + df) - 1  # note: equal to -df ~ 1e-10 if df**2 < 0.5e-10
    DT_HISTORY.clear()  # restart collection

    # Cap correction if positive value supplied
    df_applied = df_toapply
    if cap > 0:
        df_applied = max(-cap, min(cap, df_toapply))

    logger.info(
        "Correction: %5.1f ppb (%5.1f ppb @ %s)",
        df_applied*1e9, df_toapply*1e9, epochs[-1],
    )
    return df_applied

def main():
    script_name = Path(sys.argv[0]).name
    parser = configargparse.ArgumentParser(
        add_config_file_help=False,
        default_config_files=[f"{script_name}.default.conf"],
        description=parse_docstring_description(__doc__),
        formatter_class=ArgparseCustomFormatter,
        add_help=False,
    )

    # Boilerplate
    pgroup_config = parser.add_argument_group("display/configuration")
    pgroup_config.add_argument(
        "-h", "--help", action="store_true",
        help="Show this help message and exit")
    pgroup_config.add_argument(
        "-v", "--verbosity", action="count", default=0,
        help="Specify debug verbosity, e.g. -vv for more verbosity")
    pgroup_config.add_argument(
        "-L", "--logging", metavar="",
        help="Log to file, if specified. Log level follows verbosity.")
    pgroup_config.add_argument(
        "--config", metavar="file", is_config_file_arg=True,
        help="Path to configuration file")
    pgroup_config.add_argument(
        "--save", metavar="", is_write_out_config_file_arg=True,
        help="Path to configuration file for saving, then immediately exit")

    # freqcalc parameters
    pgroup = parser.add_argument_group("freqservo parameters")
    pgroup.add_argument(
        "-i", metavar="ignore", type=int, default=5,
        help="Ignore initial epochs during calculation (default: %(default)d)")
    pgroup.add_argument(
        "-a", metavar="average", type=int, default=3,
        help="Number of epochs to average servoed time over (default: %(default)d)")
    pgroup.add_argument(
        "-s", metavar="separation", type=int, default=12,
        help="Separation width for frequency calculation, epoch units (default: %(default)d)")
    pgroup.add_argument(
        "-c", metavar="cap", type=float, default=100e-9,
        help="Correction limit/cap per frequency update, disable with 0 (default: 100e-9)")

    # costream parameters
    pgroup = parser.add_argument_group("costream parameters")
    pgroup.add_argument(
        "-n", metavar="logfile5",
        help="Optional path to pipe of costream/freq outputs (default: stdin)")
    pgroup.add_argument(
        "-V", metavar="level", type=int, choices=(4,5),
        help="Output format of costream (default: '<epoch:str> <dt:float_ns>\\n')")
    pgroup.add_argument(
        "-e", action="store_true",
        help="Echo costream output back into stderr")

    # freqcd parameters
    pgroup = parser.add_argument_group("freqcd parameters")
    pgroup.add_argument(
        "-F", metavar="freqfilename",
        help="Path to freqcd frequency correction pipe (default: stdout)")
    pgroup.add_argument(
        "-f", metavar="freqcorr", type=int, default=0,
        help="freqcd initial frequency offset, used if '-u' not supplied (default: %(default)d)")
    pgroup.add_argument(
        "-d", action="store_true",
        help="freqcd in 'decimal mode' (default: False)")
    pgroup.add_argument(
        "-u", action="store_true",
        help="freqcd in 'update mode' (default: False)")

    # Parse arguments
    args = parser.parse_args()

    # Check whether options have been supplied, and print help otherwise
    if args.help:
        parser.print_help(sys.stderr)
        sys.exit(1)

    # Set logging level and log arguments
    if args.logging is not None:
        handler = logging.FileHandler(filename=args.logging, mode="w")
        handler.setFormatter(
            logging.Formatter(
                fmt=_LOGGING_FMT,
                    datefmt="%Y%m%d_%H%M%S",
                style="{",
            )
        )
        logger.addHandler(handler)

    # Set logging level
    levels = [logging.WARNING, logging.INFO, logging.DEBUG]
    verbosity = min(args.verbosity, len(levels)-1)
    logger.setLevel(levels[verbosity])
    logger.debug("%s", args)

    # Check arguments
    assert args.i > 0
    assert args.a > 0
    assert args.s > 0

    # Mainloop
    try:
        # Open relevant I/O pipes
        if args.n is not None:
            cin = open(args.n, "r")
        else:
            cin = sys.stdin
        if args.F is not None:
            cout = open(args.F, "w")
        else:
            cout = sys.stdout

        # Set initial frequency difference
        df = args.f
        if not args.d:
            df = args.f * FCORR_BTOD  # convert 2^-34 -> 0.1ppb
        df = df * 1e-10  # convert 0.1ppb -> 1

        # Restarts read if mainloop terminates, e.g. no data on read
        while True:
            line = cin.readline()
            if line == "":
                logger.debug("No new data, waiting...")
                time.sleep(0.5)
                continue

            # Derive frequency shift from costream
            epoch, dt = parse(line, format=args.V)
            if args.e:
                print(line, end="", file=sys.stderr)
            df1 = evaluate_freqshift(
                epoch, dt,
                ignore=args.i, average=args.a, separation=args.s, cap=args.c,
            )
            if df1 is None:
                logger.debug("No changes to frequency submitted.")
                continue

            # If 'update mode' is enabled on freqcd, we simply pass the correction
            # directly. Otherwise, we manually calculate the new absolute frequency.
            # Note: This is effectively the inverse of freqcd 'update mode'.
            if args.u:
                df = df1
            else:
                df = (1 + df)*(1 + df1) - 1

            # Convert and write new frequency
            _df = df * 1e10  # convert 1 -> 0.1ppb
            if not args.d:
                _df = _df * FCORR_DTOB  # convert 0.1ppb -> 2^-34
            _df = round(_df)
            logger.debug("Wrote '%d' to output", _df)
            cout.write(f"{_df}\n")
            cout.flush()  # necessary to avoid trapping df in buffer
                          # alternatively, use 'open(pipe, "wb", 0)'

    except KeyboardInterrupt:
        logger.debug("User interrupted.")

    finally:
        if args.n is not None:
            cin.close()
        if args.F is not None:
            cout.close()

if __name__ == "__main__":
    main()
