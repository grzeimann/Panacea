from __future__ import annotations

from astropy.io import fits


def test_sample_dataset_files_exist_and_shape(sample_dataset):
    # Ensure all expected files exist and have the tiny shape
    expected = ["bias", "flat", "arc", "twilight", "science"]
    for key in expected:
        path = sample_dataset[key]
        assert path.exists(), f"Missing {key} file"
        with fits.open(path) as hdul:
            data = hdul[0].data
            assert data is not None and data.ndim == 2
            ny, nx = data.shape
            # Small images only, to keep repo light and tests fast
            assert ny <= 128 and nx <= 256
            # Basic header sanity
            hdr = hdul[0].header
            assert hdr.get("INSTRUME") == "LRS2"
            assert "CHANNEL" in hdr
            assert hdr.get("EXPTIME") is not None


def test_arc_has_bright_columns(sample_dataset):
    # Arc frame should have some columns significantly brighter than median
    path = sample_dataset["arc"]
    with fits.open(path) as hdul:
        img = hdul[0].data
    import numpy as np

    col_sums = img.sum(axis=0)
    med = np.median(col_sums)
    mad = np.median(np.abs(col_sums - med)) + 1e-6
    thresh = med + 20 * mad
    assert (col_sums > thresh).any()
