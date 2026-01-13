from __future__ import annotations

from io import StringIO

import numpy as np
from astropy.table import Table

from panacea.utils import read_arc_lines, get_config_file


def test_read_arc_lines_stringio_basic():
    data = """
    # comment
    5460.735  766   0.71  Hg
    5085.822  456   1.00  Cd
    5769.598  1022  0.14  Hg
    """
    t = read_arc_lines(StringIO(data))
    assert isinstance(t, Table)
    assert list(t.colnames) == ["col1", "col2", "col3", "col4"]
    assert len(t) == 3
    # spot-check types
    assert np.isclose(t[0]["col1"], 5460.735)
    assert np.isclose(t[1]["col2"], 456)
    assert isinstance(t[2]["col4"], str)


def test_dar_table_is_three_column_numeric():
    # Use a packaged DAR file and verify we can parse 3 numeric columns per data row
    traversable = get_config_file("dar_BL.dat")
    rows = []
    with traversable.open("r") as f:
        for raw in f:
            line = raw.strip()
            if not line or line.lower().startswith("wave"):
                continue
            # skip separator-only lines
            if set(line) <= set("- "):
                continue
            parts = line.split()
            if len(parts) < 3:
                continue
            try:
                rows.append(tuple(float(x) for x in parts[:3]))
            except Exception:
                continue
    assert len(rows) > 0, "No numeric rows parsed from DAR table"
    # ensure shapes
    assert all(len(r) == 3 for r in rows)


essential_sky_files = [
    "orange_skylines.dat",
    "red_skylines.dat",
    "farred_skylines.dat",
]


def test_skylines_files_first_value_float():
    # Ensure each skylines file has at least one parseable wavelength entry
    for name in essential_sky_files:
        traversable = get_config_file(name)
        found = False
        with traversable.open("r") as f:
            for raw in f:
                line = raw.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split()
                try:
                    float(parts[0])
                    found = True
                    break
                except Exception:
                    continue
        assert found, f"No float wavelength found in {name}"
