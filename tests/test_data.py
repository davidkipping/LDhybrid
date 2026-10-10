"""Checks on downloaded data (skipped when data/ is absent)."""

import numpy as np
import pytest
from astropy.io import fits

from ldhybrid import bands
from ldhybrid.grids import phoenix
from ldhybrid.paths import RAW

pytestmark = pytest.mark.data

RAW_5800 = RAW / "phoenix" / "Z-0.0" / phoenix.file_name(5800, 4.5, "-0.0")


def test_name_roundtrip():
    n = phoenix.file_name(5800, 4.5, "-0.0")
    assert n == "lte05800-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-SpecInt.fits"
    assert phoenix.parse_name(n) == dict(teff=5800, logg=4.5, mh=-0.0, alpha=0.0)
    assert phoenix.parse_name(phoenix.file_name(3000, 5.5, "+0.5"))["mh"] == 0.5


@pytest.mark.skipif(not RAW_5800.exists(), reason="raw P2-subset file not downloaded")
def test_reduced_matches_raw_band_integration():
    s = phoenix.reduce_file(RAW_5800)
    with fits.open(RAW_5800) as h:
        I, mu = h[0].data, h["MU"].data
    lam = 500.0 + np.arange(I.shape[1])
    assert np.array_equal(s.mu, mu)
    # J top-hat on whole pixels, trapezoid of I·λ directly from the raw file
    m = (lam >= 11000) & (lam <= 14000)
    direct = np.trapezoid(I[:, m] * lam[m], lam[m], axis=1)
    red = s.band([10999.0, 11000.0, 14000.0, 14001.0], [0.0, 1.0, 1.0, 0.0])
    assert np.max(np.abs(red / red[-1] - direct / direct[-1])) < 5e-5
    # total photon integral over the binned range is conserved (float32 storage)
    e0, e1 = s.lam_edges[0], s.lam_edges[-1]
    k = (lam >= e0) & (lam <= e1)
    tot = np.trapezoid(I[:, k] * lam[k], lam[k], axis=1)
    assert np.allclose(s.P.astype(float).sum(0), tot, rtol=2e-4)


@pytest.mark.parametrize("name", bands.BANDS)
def test_band_curves_load(name):
    lam, S = bands.response(name)
    assert np.all(np.diff(lam) >= 0) and np.all(S >= 0) and S.max() > 0


MPS1 = RAW / "mpsatlas" / "set1.zip"


@pytest.mark.skipif(not MPS1.exists(), reason="MPS-ATLAS Set 1 not downloaded")
@pytest.mark.parametrize("band", ["Kepler", "TESS", "CHEOPS"])
def test_mpsatlas_band_profiles_match_kostogryz(band):
    """Our photon-weighted integration of the archived MPS-ATLAS spectra
    reproduces the published CLVs (J/A+A/666/A60 tablee3) to ~1e-3."""
    from ldhybrid import reference as ref
    from ldhybrid.grids import mpsatlas
    clv = ref.kostogryz_clv("tablee3")
    lam, S = bands.response(band)
    for mh, teff, logg in [(0.0, 5800, 4.5), (0.0, 3500, 5.0), (0.5, 4500, 4.6)]:
        s = mpsatlas.load("set1", mh, teff, logg)
        I = s.band(lam, S)
        y = np.interp(ref.KOSTOGRYZ_MU, s.mu, I / I[-1])
        k = np.flatnonzero((clv["mh"] == mh) & (clv["teff"] == teff)
                           & (clv["logg"] == logg) & (clv["band"] == band))[0]
        assert np.max(np.abs(y - clv["y"][k])) < 1.5e-3
