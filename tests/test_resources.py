from __future__ import annotations

import io

from panacea.utils import get_config_file, read_arc_lines


def test_packaged_resource_exists_and_is_readable():
    # Pick one of the known packaged files
    traversable = get_config_file('lines_uv.dat')
    assert traversable is not None
    assert hasattr(traversable, 'is_file') and traversable.is_file()

    # Open and read a few bytes to ensure packaging works in both dev and installed modes
    with traversable.open('r') as f:
        head = f.read(128)
        assert isinstance(head, str)
        assert len(head) > 0


def test_read_arc_lines_parses_minimal_valid_content():
    # Create a small in-memory file compatible with read_arc_lines expectations
    content = """
# comment line
4000.0  10.0  0.5  Hg
5000.0  20.0  1.0  Ne
bad line
6000    30  0.2  Ar
"""
    tbl = read_arc_lines(io.StringIO(content))
    # Should parse 3 valid data rows, 4 named columns
    assert tbl.shape == (3, 4)
    assert list(tbl.colnames) == ['col1', 'col2', 'col3', 'col4']
    # Check a value
    assert tbl['col1'][0] == 4000.0
