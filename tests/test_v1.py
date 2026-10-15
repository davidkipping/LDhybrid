"""The v1.0 law definitions are frozen (Gate G2): any change here is a new version."""

import numpy as np

from ldhybrid import laws, v1


def test_frozen_poles():
    assert v1.HYBRID2_EPS == 0.36
    assert laws.get("hybrid2").shapes[1].eps.tolist() == [0.36]
    eps4 = sorted({float(e) for s in laws.get("hybrid4").shapes for e in s.eps})
    assert np.allclose(eps4, [0.02860121, 0.38582418], rtol=1e-6)
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


def _g_derivs(L, w, xi):
    """dI/dxi, d2I/dxi2, d3I/dxi3 (xi = mu^2) by central differences of the
    exact slope."""
    h = 1e-5
    g = -L.dT(xi) @ w
    gp = -(L.dT(xi + h) - L.dT(xi - h)) @ w / (2 * h)
    gpp = -(L.dT(xi + h) - 2 * L.dT(xi) + L.dT(xi - h)) @ w / h ** 2
    return g, gp, gpp


def test_hybrid4_basis_is_the_curvature_simplex():
    L = laws.get("hybrid4")
    xi = np.linspace(1e-4, 1 - 1e-4, 801)
    rng = np.random.default_rng(1)
    # every point of the simplex: positive, monotonic, concave in mu^2, convex slope
    W = rng.dirichlet(np.ones(5), 3000)[:, :4]
    for w in W:
        g, gp, gpp = _g_derivs(L, w, xi)
        assert 1 - w.sum() >= -1e-12 and g.min() >= -1e-9 and gp.max() <= 1e-6 and gpp.min() >= -1e-3
    # every vertex sits on the boundary: S3, S4 have zero limb intensity and
    # zero slope and curvature at the disc centre
    for k in (2, 3):
        w = np.eye(4)[k]
        g, gp, _ = _g_derivs(L, w, np.array([1.0]))         # analytic beyond 1
        assert abs(1 - w.sum()) < 1e-12 and abs(g[0]) < 1e-6 and abs(gp[0]) < 1e-4
    # just outside each facet w_k = 0 a condition fails
    for k in range(4):
        w = np.full(4, 0.1)
        w[k] = -0.01
        g, gp, gpp = _g_derivs(L, w, xi)
        assert g.min() < 0 or gp.max() > 0 or gpp.min() < 0
