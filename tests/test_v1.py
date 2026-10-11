"""The v1.0 law definitions are frozen (Gate G2): any change here is a new version."""

import numpy as np

from ldhybrid import laws, v1


def test_frozen_poles():
    assert v1.HYBRID2_EPS == 0.36
    assert laws.get("hybrid2").shapes[1].eps.tolist() == [0.36]
    assert np.allclose(laws.get("hybrid4").shapes[2].eps.tolist() + laws.get("hybrid4").shapes[3].eps.tolist(),
                       [0.02860121, 0.38582418], rtol=1e-6)
    eps5 = [s.eps[0] for s in laws.get("hybrid5").shapes[2:]]
    assert np.allclose(eps5, [0.00159655, 0.04127504, 0.38743472], rtol=1e-6)


def test_hybrid2_triangle_vertices_are_the_physical_boundary():
    import squishierplanet.laws as sl
    L = laws.get("hybrid2")
    Vc, Vl = sl.hybrid2_vertices(v1.HYBRID2_EPS)
    mu = np.linspace(0, 1, 2001)
    for V in (Vc, Vl):
        y = L.intensity(V, mu)
        assert abs(y[0]) < 1e-12                       # limb intensity zero
        assert np.all(np.diff(y) >= -1e-12)            # monotone
    assert np.allclose(Vc, [-0.1246, 1.1246], atol=1e-4)
    assert np.allclose(Vl, [1.2010, -0.2010], atol=1e-4)


def test_method_constants():
    assert (v1.MEASURE, v1.RESAMPLE, v1.FCM, v1.LIMB) == ("r", "exact", False, "M1")
