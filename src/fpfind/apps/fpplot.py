#!/usr/bin/env python3
"""NOT FOR PRODUCTION USE - but in main branch for ease of maintenance ;)

Changelog:
    2024-02-02, Justin: Init
"""

import argparse
import logging
import sys

import kochen  # v0.2024.4
import kochen.logging
import kochen.scriptutil
import matplotlib.pyplot as plt
import numpy as np
from S15lib.g2lib.g2lib import histogram

from fpfind.lib.parse_timestamps import read_a1_overlapping
from fpfind.lib.utils import (
    get_first_overlapping_epoch,
    get_timestamp_pattern,
    normalize_timestamps,
)

_ENABLE_BREAKPOINT = False
logger = logging.getLogger(__name__)


def plotter(
    alice,
    bob,
    freq,
    time,
    width,
    resolution,
    window=800,
    save=False,
    normalize=False,
    decimation=1,
):
    bob = (bob - time) / (1 + freq * 1e-6)
    ys, xs, stats = histogram(
        alice,
        bob,
        duration=width,
        resolution=resolution,
        statistics=True,
        window=window,
    )
    # Custom breakpoint for experimentation
    if _ENABLE_BREAKPOINT:
        globals().update(locals())  # write all local variables to global scope
        raise  # noqa: PLE0704

    if decimation > 1:
        ys = ys[::decimation]
        xs = xs[::decimation]

    elapsed = (alice[-1] - alice[0]) * 1e-9
    _center = ys[500:1500]
    _sides = np.array(list(ys[:500]) + list(ys[1500:]))
    coincidences = np.sum(_center) - np.sum(_sides) * len(_center) / len(_sides)

    s1 = len(alice) / elapsed
    s2 = len(bob) / elapsed
    c = coincidences / elapsed
    print("Totals:")
    print(f"  S1: {len(alice):d}")
    print(f"  S2: {len(bob):d}")
    print(f"  C:  {coincidences:.0f}")
    print(f"Time elapsed: {elapsed:.3f}s")
    print(f"  s1: {s1:.0f}")
    print(f"  s2: {s2:.0f}")
    print(f"  sg: {(len(alice)*len(bob))**0.5/elapsed:.0f}")
    print(f"  c:  {c:.0f}")
    print(f"Efficiency: {100*c/np.sqrt(s1*s2):.1f}%")

    # Normalize?
    if normalize:
        accs = s1 * s2 * elapsed * resolution * 1e-9
        yerrs = np.sqrt(ys) / accs
        ys = ys / accs
        print(f"Normalizing by background = {accs}")
        plt.errorbar(xs, ys, yerrs, fmt=".", markersize=3, elinewidth=1)

    else:
        plt.plot(xs, ys, linewidth=1)

    plt.xlabel("Delay, $\\tau$ (ns)")
    plt.ylabel("$g^{(2)}(\\tau)$")
    plt.tight_layout()
    if save:
        plt.xlim(np.min(xs), np.max(xs))
        plt.tight_layout()
        plt.savefig(save)
    plt.show()
    globals().update(locals())
    raise  # noqa: PLE0704


def attenuate(ts, transmission=1):
    if transmission >= 1:
        return ts
    mask = np.random.random(size=ts.size) < transmission
    return ts[mask]


def main():  # noqa: PLR0915
    global _ENABLE_BREAKPOINT  # noqa: PLW0603
    parser = kochen.scriptutil.generate_default_parser(__doc__, "fpfind")

    # Boilerplate
    # fmt: off
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

    # Epoch importing arguments
    pgroup_chselect = parser.add_argument_group("channel selection")
    pgroup_chselect.add_argument(
        "-m", "--reference-pattern", metavar="", type=int,
        help="Pattern mask for selecting detector events from low-count side")
    pgroup_chselect.add_argument(
        "-M", "--target-pattern", metavar="", type=int,
        help="Pattern mask for selecting detector events from high-count side")

    # Plotting parameters
    pgroup = parser.add_argument_group("plotting")
    pgroup.add_argument(
        "--df", type=float, default=0.0,
        help="Specify clock skew, in units of ppm (default: %(default)f)")
    pgroup.add_argument(
        "--dt", type=float, default=0.0,
        help="Specify time delay, in units of ns (default: %(default)f)")
    pgroup.add_argument(
        "--width", type=float, default=1000,
        help="Specify one-sided width of histogram, in units of ns (default: %(default)f)")  # noqa: E501
    pgroup.add_argument(
        "-r", "--resolution", "--final-res", type=float, default=1,
        help="Specify resolution of histogram, in units of ns (default: %(default)f)")
    pgroup.add_argument(
        "--duration", type=float, default=1,
        help="Specify duration of timestamps to use, units of s (default: %(default)f)")
    pgroup.add_argument(
        "--save-plot",
        help="Specify filename to save the plot to")
    pgroup.add_argument(
        "--normalize", action="store_true",
        help="Normalize g(2) plot")
    pgroup.add_argument(
        "--decimation", type=int, default=1,
        help="Take only every N-th point, for clarity (default: %(default)d)")
    # fmt: on

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
            args.sendfiles,
            args.t1files,
            first_epoch=args.first_epoch,
            return_length=True,
        )
        if available_epochs < args.num_epochs + args.skip_epochs:
            logger.warning("  Insufficient epochs")

        # Read epochs
        alice, aps = get_timestamp_pattern(
            args.sendfiles, "T2", first_epoch, args.skip_epochs, args.num_epochs
        )
        bob, bps = get_timestamp_pattern(
            args.t1files, "T1", first_epoch, args.skip_epochs, args.num_epochs
        )

    elif args.target is not None and args.reference is not None:
        logger.info("  Reading from timestamp files...")
        _is_reading_ts = True
        (alice, bob), (aps, bps) = read_a1_overlapping(
            args.reference,
            args.target,
            legacy=args.legacy,
            duration=args.duration + args.skip_duration,
        )

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

    alice = attenuate(alice, 1)
    bob = attenuate(bob, 1)

    plotter(
        alice,
        bob,
        args.df,
        args.dt,
        2 * args.width,
        resolution=args.resolution,
        save=args.save_plot,
        normalize=args.normalize,
        decimation=args.decimation,
    )


if __name__ == "__main__":
    main()
