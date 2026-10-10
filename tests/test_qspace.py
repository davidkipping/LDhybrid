import numpy as np
import squishierplanet.laws as sl

from ldhybrid import qspace, v1


def test_hybrid2_roundtrip():
    q = np.random.default_rng(0).uniform(size=(200, 2))
    w = sl.hybrid2_from_q(q[:, 0], q[:, 1], v1.HYBRID2_EPS)
    back = np.array([qspace.hybrid2_to_q(x) for x in w])
    assert np.allclose(back, q, atol=1e-10)


def test_simplex_roundtrip():
    for n in (4, 5):
        q = np.random.default_rng(n).uniform(size=(200, n))
        w = sl.simplex_from_q(q)
        back = np.array([qspace.simplex_to_q(x) for x in w])
        assert np.allclose(back, q, atol=1e-9)
