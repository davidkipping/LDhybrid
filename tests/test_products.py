"""P3 product checks on one model (skipped without data)."""

import numpy as np
import pytest

from ldhybrid import bands, qspace, spectra
from ldhybrid.paths import REDUCED

pytestmark = pytest.mark.data
HAVE = (REDUCED / "phoenix" / "Z-0.0").exists()


@pytest.mark.skipif(not HAVE, reason="PHOENIX not reduced")
def test_model_products_physical_and_complete():
    from ldhybrid.grids import phoenix
    from ldhybrid.products import LAWS, model_products
    bl = ["TESS", "HST_WFC3_G141", "JWST_NIRSpec_G395H"]
    r = model_products(phoenix.load(5800, 4.5, "-0.0"), bl)
    assert r["covered"].all() and np.all((r["r_c"] > 0.99) & (r["r_c"] < 1))
    for law in ("hybrid2", "hybrid4", "hybrid5"):
        assert r[f"{law}|r|physical"].all()
        assert np.all(r[f"{law}|r|lc_ppm"] < {"hybrid2": 60, "hybrid4": 15, "hybrid5": 10}[law])
        for w in r[f"{law}|r|w"]:
            q = qspace.to_q(law, w)
            assert np.all((q >= -1e-9) & (q <= 1 + 1e-9))
    assert set(LAWS) >= {"hybrid2", "hybrid4", "hybrid5"}


@pytest.mark.skipif(not HAVE, reason="PHOENIX not reduced")
def test_coefficient_spectrum_closure():
    from ldhybrid.grids import phoenix
    s = phoenix.load(4600, 4.5, "-0.0")
    cs = spectra.model_spectrum(s)
    for b in ("TESS", "JWST_NIRISS_SOSS_o1"):
        lam, S = bands.response(b)
        for law in spectra.SPECTRA_LAWS:
            w, _ = cs.band_coefficients(lam, S, law)
            ref = spectra.band_fit_same_settings(s, lam, S, law, cs.r_c)
            assert np.allclose(w, ref, atol=1e-12)
    # rebinning is exact for a top-hat aligned with the merged bins
    e = cs.lam_edges
    top = ([e[300] - 1e-6, e[300], e[900], e[900] + 1e-6], [0, 1, 1, 0])
    w1, _ = cs.band_coefficients(*top, "hybrid5")
    w3, _ = cs.rebin(3).band_coefficients(*top, "hybrid5")
    assert np.allclose(w1, w3, atol=1e-12)
