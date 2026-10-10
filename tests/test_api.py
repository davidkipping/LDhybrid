import numpy as np
import pytest

import ldhybrid as ld
from ldhybrid import bands
from ldhybrid.paths import DATA

pytestmark = pytest.mark.data
HAVE = (DATA / "products" / "mps2.h5").exists() and (DATA / "products" / "spectra").exists()


@pytest.mark.skipif(not HAVE, reason="products not built")
def test_coefficients_and_correction():
    r = ld.coefficients("hybrid5", 5600, 4.4, 0.1, "TESS", grid="mps2")
    assert r["physical"] and np.all((r["q"] >= 0) & (r["q"] <= 1))
    rc = ld.coefficients("hybrid5", 5600, 4.4, 0.1, "TESS", grid="mps2", correction="B100G")
    assert not np.allclose(r["w"], rc["w"]) and rc["physical"]
    law = ld.to_law("hybrid5", r["w"])
    assert np.isclose(law.intensity(np.array([1.0]))[0], 1.0)
    with pytest.raises(ValueError):
        ld.coefficients("hybrid5", 9000, 4.4, 0.0, "TESS", grid="mps2")


@pytest.mark.skipif(not HAVE, reason="products not built")
def test_custom_band_close_to_table():
    lam, S = bands.response("HST_WFC3_G141")
    c = ld.coefficients_custom("hybrid4", 5800, 4.5, 0.0, lam, S)
    t = ld.coefficients("hybrid4", 5800, 4.5, 0.0, "HST_WFC3_G141", grid="phoenix")
    # one r_c per model vs per band: the weakly determined inner-pole weight can
    # differ by a few 0.01, so compare light curves (P3: 2-10 ppm)
    from ldhybrid import laws
    from ldhybrid import lightcurve as lc
    L = laws.get("hybrid4")
    d = lc.lc_errors(lambda m: L.intensity(c["w"], m), lambda m: L.intensity(t["w"], m))["max_ppm"]
    assert c["physical"] and d < 15
