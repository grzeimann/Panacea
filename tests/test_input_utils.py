from __future__ import annotations

import argparse as ap
import types
import pytest

from panacea import input_utils as iu


def test_setup_parser_defaults_and_flags():
    parser = iu.setup_parser()
    # Check that key options exist and defaults are as expected
    ns = parser.parse_args([])
    assert ns.start_date is None
    assert ns.end_date is None
    assert ns.date_length is None
    assert ns.rootdir == '/work/03946/hetdex/maverick'
    assert ns.instrument == 'virus'

    # Parse with a mixture of options
    ns = parser.parse_args(["--start_date", "20250101", "--end_date", "20250103"])  # 2-day open interval
    # Add a logger required by set_daterange
    ns.log = iu.setup_logging("test")
    ns = iu.set_daterange(ns)
    # Expect dates 2025-01-01 and 2025-01-02 (end exclusive)
    assert len(ns.daterange) == 2
    assert str(ns.daterange[0]) == "2025-01-01"
    assert str(ns.daterange[1]) == "2025-01-02"


def test_setup_parser_length_forward_and_backward():
    parser = iu.setup_parser()

    # Forward range from start_date for N days
    ns = parser.parse_args(["--start_date", "20240227", "--date_length", "3"])  # 27,28,29
    ns.log = iu.setup_logging("test")
    ns = iu.set_daterange(ns)
    assert [str(d) for d in ns.daterange] == ["2024-02-27", "2024-02-28", "2024-02-29"]

    # Backward range from end_date for N days (descending in function implementation)
    ns = parser.parse_args(["--end_date", "20240227", "--date_length", "2"])  # 27, 26
    ns.log = iu.setup_logging("test")
    ns = iu.set_daterange(ns)
    # Order is end_date - 0 days, then -1 day according to implementation
    assert [str(d) for d in ns.daterange] == ["2024-02-27", "2024-02-26"]


def test_setup_basic_parser_defaults_and_parse():
    parser = iu.setup_basic_parser()
    ns = parser.parse_args([])
    assert ns.date is None
    assert ns.observation is None
    assert ns.exposure_number is None
    assert ns.rootdir == '/work/03946/hetdex/maverick'
    assert ns.instrument == 'lrs2'
    assert ns.ifuslot == '066'
    assert ns.side == 'L'

    # Provide a full set of arguments
    ns = parser.parse_args([
        "-d", "20181108",
        "-o", "7",
        "-e", "10",
        "-r", "/tmp/root",
        "-in", "lrs2",
        "-i", "070",
        "-s", "R",
    ])
    assert ns.date == "20181108"
    assert ns.observation == "7"
    assert ns.exposure_number == 10
    assert ns.rootdir == "/tmp/root"
    assert ns.instrument == "lrs2"
    assert ns.ifuslot == "070"
    assert ns.side == "R"
