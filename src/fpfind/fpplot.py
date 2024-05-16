#!/usr/bin/env python3
"""NOT FOR PRODUCTION USE - but in main branch for ease of maintenance ;)

Changelog:
    2024-02-02, Justin: Init
"""

import logging
import sys

import argparse
import matplotlib.pyplot as plt

import kochen.scriptutil
import kochen.logging
from S15lib.g2lib.g2lib import histogram

from fpfind.lib.parse_timestamps import read_a1
from fpfind.lib.utils import (
    get_first_overlapping_epoch, get_timestamp,
    normalize_timestamps, slice_timestamps,
)

from fpfind.lib.utils import normalize_timestamps
from fpfind.lib.parse_timestamps import read_a1_overlapping

_ENABLE_BREAKPOINT = False
logger = logging.getLogger(__name__)

def plotter(alice, bob, freq, time, width, resolution, window=800, save=False):
    bob = (bob - time) / (1 + freq*1e-6)
    ys, xs, stats = histogram(
        alice, bob, duration=width, resolution=resolution,
        statistics=True, window=window,
    )
    # Custom breakpoint for experimentation
    if _ENABLE_BREAKPOINT:
        globals().update(locals())  # write all local variables to global scope
        raise

    elapsed = (alice[-1] - alice[0]) * 1e-9
    coincidences = stats.signal.sum() - (stats.background.sum() * stats.signal.size / stats.background.size)

    print(f"Totals:")
    print(f"  S1: {len(alice):d}")
    print(f"  S2: {len(bob):d}")
    print(f"  C:  {coincidences:.0f}")
    print(f"Time elapsed: {elapsed:.3f}s")
    print(f"  s1: {len(alice)/elapsed:.0f}")
    print(f"  s2: {len(bob)/elapsed:.0f}")
    print(f"  sg: {(len(alice)*len(bob))**0.5/elapsed:.0f}")
    print(f"  c:  {coincidences/elapsed:.0f}")

    plt.plot(xs, ys, linewidth=1)
    plt.xlabel("Delay (ns)")
    plt.ylabel("g(2)")
    plt.tight_layout()
    if save:
        plt.savefig(save)
    plt.show()


def main():
    global _ENABLE_BREAKPOINT
    parser = kochen.scriptutil.generate_default_parser(__doc__)

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
        "--quiet", action="store_true",
        help="Suppress errors, but will not block logging")
    pgroup_config.add_argument(
        "--config", metavar="", is_config_file_arg=True,
        help="Path to configuration file")
    pgroup_config.add_argument(
        "--save", metavar="", is_write_out_config_file_arg=True,
        help="Path to configuration file for saving, then immediately exit")
    pgroup_config.add_argument(
        "--experiment", action="store_true",
        help=argparse.SUPPRESS)

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
        help="Specify initial duration to skip, in seconds (default: %(default)s)")

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
        help="Specify filename of first overlapping epoch, optional")
    pgroup_ep.add_argument(
        "-n", "--num-epochs", metavar="", type=int, default=1,
        help="Specify number of epochs to import (default: %(default)d)")
    pgroup_ep.add_argument(
        "-z", "--skip-epochs", metavar="", type=int, default=0,
        help="Specify number of initial epochs to skip (default: %(default)d)")

    # Plotting parameters
    pgroup = parser.add_argument_group("plotting")
    pgroup.add_argument(
        "--freq", type=float, default=0.0,
        help="Specify clock skew, in units of ppm (default: %(default)f)")
    pgroup.add_argument(
        "--time", type=float, default=0.0,
        help="Specify time delay, in units of ns (default: %(default)f)")
    pgroup.add_argument(
        "--width", type=float, default=1000,
        help="Specify width of histogram, in units of ns (default: %(default)f)")
    pgroup.add_argument(
        "--resolution", type=float, default=1,
        help="Specify resolution of histogram, in units of ns (default: %(default)f)")
    pgroup.add_argument(
        "--duration", type=float,
        help="Specify duration of timestamps to use, units of s")
    pgroup.add_argument(
        "--save-plot",
        help="Specify filename to save the plot to")


    # Parse arguments and configure logging
    args = kochen.scriptutil.parse_args_or_help(parser, ignore_unknown=True)
    kochen.logging.set_default_handlers(logger, file=args.logging)
    kochen.logging.set_logging_level(logger, args.verbosity)
    logger.debug("%s", args)

    # Set experimental mode
    if args.experiment:
        _ENABLE_BREAKPOINT = True

    # Obtain timestamps needed for fpplot
    #   alice: low count side - chopper - HeadT2 - sendfiles (reference)
    #   bob: high count side - chopper2 - HeadT1 - t1files
    if args.sendfiles is not None and args.t1files is not None:
        logger.info("  Reading from epoch directories...")
        _is_reading_ts = False

        # Automatically choose first overlapping epoch if not supplied manually
        first_epoch, available_epochs = get_first_overlapping_epoch(
            args.sendfiles, args.t1files,
            first_epoch=args.first_epoch, return_length=True,
        )
        if available_epochs < args.num_epochs + args.skip_epochs:
            logger.warning("  Insufficient epochs")

        # Read epochs
        alice = get_timestamp(
            args.sendfiles, "T2",
            first_epoch, args.skip_epochs, args.num_epochs)
        bob = get_timestamp(
            args.t1files, "T1",
            first_epoch, args.skip_epochs, args.num_epochs)

    elif args.target is not None and args.reference is not None:
        logger.info("  Reading from timestamp files...")
        _is_reading_ts = True
        (alice, bob), _ = read_a1_overlapping(
            args.reference, args.target, legacy=args.legacy,
            duration=args.duration,
        )

    else:
        logger.error("Timestamp files/epochs must be supplied with -tT/-dD")
        sys.exit(1)

    # Normalize timestamps to common time reference near start, so that
    # frequency compensation will not shift the timing difference too far
    skip = args.skip_duration if _is_reading_ts else 0
    alice, bob = normalize_timestamps(alice, bob, skip=skip)

    plotter(alice, bob, args.freq, args.time, args.width, resolution=args.resolution, save=args.save_plot)


if __name__ == "__main__":
    main()
