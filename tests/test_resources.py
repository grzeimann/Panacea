from __future__ import annotations

from importlib import resources
from pathlib import Path

import panacea
from panacea.utils import get_config_file


REQUIRED_FILES = [
    # Line lists
    "lines_uv.dat",
    "lines_orange.dat",
    # DAR table(s)
    "dar_BL.dat",
    # FPlane for astrometry
    "fplane.txt",
]


def test_lrs2config_dir_exists_and_populated():
    """Verify the packaged lrs2_config directory exists and has content."""
    lrs2config = resources.files("panacea") / "lrs2_config"
    # resources.files returns a Traversable; ensure it exists and is a directory-like object
    assert lrs2config is not None

    # Ensure there is at least one file inside (any extension)
    # Convert to list via iterdir() if supported
    children = list(lrs2config.iterdir())
    assert len(children) > 0, "lrs2_config appears to be empty"


def test_required_resource_files_present():
    """Check that a few key resource files can be located via get_config_file."""
    for name in REQUIRED_FILES:
        traversable = get_config_file(name)
        assert traversable is not None, f"get_config_file returned None for {name}"
        # Traversable should represent an existing file packaged with the module
        assert traversable.name == name
        # Some Traversable backends do not have .is_file(); try opening instead
        with traversable.open("rb") as fh:
            blob = fh.read(64)
            assert isinstance(blob, (bytes, bytearray))
