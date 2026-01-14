# -*- coding: utf-8 -*-
"""
Created on Wed Nov 29 09:38:06 2017

@author: gregz
"""

import sys
import datetime
import logging

import argparse as ap
from datetime import datetime as dt


def setup_parser():
    """Build a CLI parser for selecting a date range and instrument paths.

    This parser is intended for utilities that iterate over multiple dates.
    It supports specifying an explicit start and end date, or a start/end
    date plus a day-count span to construct a range. Paths and instrument
    name can also be provided.

    Returns:
        argparse.ArgumentParser: Configured parser with the following options:
        - --start_date YYYYMMDD
        - --end_date YYYYMMDD
        - --date_length N (days)
        - --rootdir PATH
        - --instrument NAME
    """
    parser = ap.ArgumentParser(add_help=True)

    parser.add_argument(
        "-sd", "--start_date", help="""Start date YYYYMMDD""", type=str, default=None
    )

    parser.add_argument(
        "-ed", "--end_date", help="""End date YYYYMMDD""", type=str, default=None
    )

    parser.add_argument(
        "-dl",
        "--date_length",
        help="""Number of days to include""",
        type=int,
        default=None,
    )

    parser.add_argument(
        "-r",
        "--rootdir",
        help="""Root directory for date tree""",
        type=str,
        default="/work/03946/hetdex/maverick",
    )

    parser.add_argument(
        "-in",
        "--instrument",
        help="""Instrument name (e.g., virus)""",
        type=str,
        default="virus",
    )

    return parser


def setup_basic_parser():
    """Build a simple CLI parser for a single date/observation selection.

    This parser is intended for commands that operate on one observation or
    exposure at a time. It allows specifying the observation date, the
    observation number and exposure number, along with root directory, the
    instrument, IFU slot, and channel side.

    Returns:
        argparse.ArgumentParser: Configured parser with options including
        --date, --observation, --exposure_number, --rootdir, --instrument,
        --ifuslot, and --side.
    """
    parser = ap.ArgumentParser(add_help=True)

    parser.add_argument(
        "-d", "--date", help="""Observation date YYYYMMDD""", type=str, default=None
    )

    parser.add_argument(
        "-o",
        "--observation",
        help="""Observation number (e.g., 7 or 0000007)""",
        type=str,
        default=None,
    )

    parser.add_argument(
        "-e",
        "--exposure_number",
        help="""Exposure number (e.g., 10)""",
        type=int,
        default=None,
    )

    parser.add_argument(
        "-r",
        "--rootdir",
        help="""Root directory for reductions""",
        type=str,
        default="/work/03946/hetdex/maverick",
    )

    parser.add_argument(
        "-in",
        "--instrument",
        help="""Instrument name (e.g., lrs2)""",
        type=str,
        default="lrs2",
    )

    parser.add_argument(
        "-i", "--ifuslot", help="""IFU slot (e.g., 066)""", type=str, default="066"
    )

    parser.add_argument(
        "-s", "--side", help="""Instrument side (e.g., L)""", type=str, default="L"
    )

    return parser


def setup_logging(logname="input_utils"):
    """Create and configure a module logger writing to stdout.

    The logger is created only once per process (subsequent calls reuse the
    existing handler). Messages are formatted with level and timestamp. The
    stream handler level is set to INFO and the logger level to DEBUG so that
    downstream code can adjust verbosity by changing the logger level.

    Args:
        logname: Name of the logger to create or retrieve. Defaults to
            'input_utils'.

    Returns:
        logging.Logger: Configured logger instance.
    """
    log = logging.getLogger("input_utils")
    if not len(log.handlers):
        fmt = "[%(levelname)s - %(asctime)s] %(message)s"
        fmt = logging.Formatter(fmt)

        level = logging.INFO

        handler = logging.StreamHandler()
        handler.setFormatter(fmt)
        handler.setLevel(level)

        log = logging.getLogger("input_utils")
        log.setLevel(logging.DEBUG)
        log.addHandler(handler)
    return log


def set_daterange(args):
    """Derive a list of dates from provided start/end and/or length options.

    Expects an argparse-like namespace containing some combination of
    start_date, end_date, and date_length. If date_length is None, both
    start_date and end_date must be provided and an inclusive range is built
    from start (inclusive) to end (exclusive). If date_length is provided, the
    function builds a forward range from start_date for N days or a backward
    range from end_date for N days. The computed list is stored on the args
    object as ``args.daterange``.

    Args:
        args: Namespace with attributes start_date (str YYYYMMDD), end_date
            (str YYYYMMDD), date_length (int or None), and log (logger with
            .error/.warning methods).

    Returns:
        The same args namespace with an added ``daterange`` attribute, a list of
        datetime.date objects representing the selected dates.
    """
    dateatt = ["start_date", "end_date"]
    if args.date_length is None:
        if args.start_date is None:
            args.log.error(
                "You must include two of the following: "
                '"start_date", "end_date", or "date_length"'
            )
            sys.exit(1)
        if args.end_date is None:
            args.log.error(
                "You must include two of the following: "
                '"start_date", "end_date", or "date_length"'
            )
            sys.exit(1)
        dates = {}
        for da in dateatt:
            dates[da] = dt(
                int(getattr(args, da)[:4]),
                int(getattr(args, da)[4:6]),
                int(getattr(args, da)[6:]),
            )

        args.daterange = [
            datetime.date.fromordinal(i)
            for i in range(dates[dateatt[0]].toordinal(), dates[dateatt[1]].toordinal())
        ]
    else:
        if args.start_date is not None and args.end_date is not None:
            args.log.warning(
                'Using "start_date" and "date_length", '
                'however, you specified "end_date" as well '
                "which will not be used."
            )
            args.end_date = None
        if args.start_date is not None:
            base = datetime.date(
                int(args.start_date[:4]),
                int(args.start_date[4:6]),
                int(args.start_date[6:]),
            )
            args.daterange = [
                base + datetime.timedelta(days=x) for x in range(0, args.date_length)
            ]

        if args.end_date is not None:
            base = datetime.date(
                int(args.end_date[:4]), int(args.end_date[4:6]), int(args.end_date[6:])
            )
            args.daterange = [
                base - datetime.timedelta(days=x) for x in range(0, args.date_length)
            ]

    return args
