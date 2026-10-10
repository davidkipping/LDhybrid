import numpy as np
import pytest

from ldhybrid.specint import SpecInt, bin_photons, log_edges


def test_log_edges_constant_resolution_and_bounds():
    e = log_edges(1000.0, 54990.0, 3000.0)
    assert e[0] == 1000.0 and e[-1] <= 54990.0
    assert np.allclose(e[1:] / e[:-1], np.exp(1 / 3000.0))


def test_bin_photons_exact_for_piecewise_linear_photon_density():
    # I·λ piecewise linear in λ -> bin integrals are exact
    lam = np.linspace(1000.0, 2000.0, 1001)
    mu = np.array([0.5, 1.0])
    f = np.vstack([1 + 1e-3 * lam, 2 - 5e-4 * lam])          # target I·λ
    I = f / lam
    edges = np.array([1000.0, 1234.5, 1500.25, 1999.9])
    P = bin_photons(lam, I, edges)
    exact = np.array([[(1 * (b - a) + 0.5e-3 * (b * b - a * a)),
                       (2 * (b - a) - 0.25e-3 * (b * b - a * a))]
                      for a, b in zip(edges[:-1], edges[1:])])
    assert np.allclose(P, exact, rtol=1e-12)


def test_bin_photons_rejects_edges_outside_samples():
    lam = np.linspace(1000.0, 2000.0, 11)
    with pytest.raises(ValueError):
        bin_photons(lam, np.ones((1, 11)), np.array([900.0, 1500.0]))


def test_band_tophat_and_coverage():
    edges = log_edges(1000.0, 5000.0, 100.0)
    s = SpecInt(edges, np.array([0.1, 1.0]), np.ones((len(edges) - 1, 2)) *
                np.diff(edges)[:, None])                     # P = bin width: flat
    lam, S = [1999.0, 2000.0, 3000.0, 3001.0], [0.0, 1.0, 1.0, 0.0]
    assert s.coverage(lam, S) == pytest.approx(1.0)
    # flat photon density: band value = ∫ S dλ
    assert s.band(lam, S) == pytest.approx(np.trapezoid(S, lam), rel=1e-12)
    with pytest.raises(ValueError):
        s.band([4000.0, 4001.0, 6000.0, 6001.0], [0.0, 1.0, 1.0, 0.0])


def test_roundtrip(tmp_path):
    e = log_edges(1000.0, 1100.0, 50.0)
    s = SpecInt(e, np.array([0.2, 1.0]), np.arange(2.0 * (len(e) - 1)).reshape(-1, 2),
                dict(teff=5000.0))
    s.save(tmp_path / "x.npz")
    t = SpecInt.load(tmp_path / "x.npz")
    assert np.array_equal(t.P, s.P) and t.meta == s.meta
