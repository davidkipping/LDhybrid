import numpy as np
import pytest

from ldhybrid import laws, predict
from ldhybrid.paths import DATA

pytestmark = pytest.mark.data
HAVE = (DATA / "products" / "phoenix.h5").exists()


def test_spam_returns_own_coefficients_for_a_quadratic_star():
    Lq = laws.quadratic()
    u = np.array([0.45, 0.18])
    for p, b in ((0.05, 0.1), (0.12, 0.7)):
        assert np.allclose(predict.spam(lambda m: Lq.intensity(u, m), p, b), u, atol=1e-8)


def test_h_prime_of_linear_law():
    L = laws.linear()
    h1, h2 = predict.h_prime(lambda m: L.intensity([0.6], m))
    assert np.isclose(h1, 1 - 0.6 / 3) and np.isclose(h2, 0.6 / 3)


@pytest.mark.skipif(not HAVE, reason="products not built")
def test_interpolation_exact_at_nodes_and_physical_between():
    import h5py
    from ldhybrid.fit import is_physical
    h = h5py.File(DATA / "products" / "phoenix.h5")
    bands = [b.decode() for b in h["bands"][()]]
    i = np.flatnonzero((h["teff"][()] == 5000) & (h["logg"][()] == 4.5) & (h["mh"][()] == 0))[0]
    w = predict.weights("phoenix", "hybrid5", "TESS", 5000, 4.5, 0.0)
    assert np.allclose(w, h["hybrid5|r|w"][i, bands.index("TESS")])
    for t in (4321, 5555, 6789):
        assert is_physical(laws.get("hybrid5"), predict.weights("phoenix", "hybrid5", "TESS", t, 4.37, -0.13))
