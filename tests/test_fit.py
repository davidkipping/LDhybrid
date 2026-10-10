import numpy as np
import pytest
import squishierplanet as sp

from ldhybrid import fit, laws
from ldhybrid import lightcurve as lc
from ldhybrid.profile import Profile, limb_cut, quadrature

MU = np.concatenate([np.linspace(0.0, 0.1, 40)[:-1], np.linspace(0.1, 1.0, 60)])
W5 = np.array([0.09, 0.16, 0.17, 0.09, 0.20])


def law_profile(name, w, mu=MU):
    L = laws.get(name)
    return Profile(mu, L.intensity(w, mu))


@pytest.mark.parametrize("measure", ["mu", "r", "area", "nodes"])
@pytest.mark.parametrize("resample", [None, 2000])
def test_recovers_exact_law(measure, resample):
    pr = law_profile("hybrid5", W5)
    for con in ("none", "simplex", "physical"):
        f = fit.fit_linear(laws.get("hybrid5"), pr, measure, resample, con)
        # resampling interpolates the profile by PCHIP: ~1e-5 near the limb
        tol = 1e-9 if resample is None or measure == "nodes" else 2e-5
        assert np.allclose(f.w, W5, atol=tol), (con, f.w)
        assert f.extra["physical"]


def test_fcm_equality_holds_and_constraints_respected():
    # a profile outside every law: a power-2 star
    pr = Profile(MU, laws.power2(0.7, 0.5, MU))
    for name in ("hybrid2", "hybrid4", "hybrid5", "claret4", "quadratic"):
        f = fit.fit_linear(laws.get(name), pr, "r", 2000, "physical", fcm=True)
        assert abs(f.flux_err) < 1e-9
        assert f.extra["physical"]


def test_physical_constraint_binds_on_unphysical_target():
    # brightening toward the limb near μ ~ 0.3: unconstrained fit is unphysical
    y = 1.0 - 0.5 * (1.0 - MU) + 0.15 * np.exp(-((MU - 0.3) / 0.05) ** 2)
    pr = Profile(MU, y)
    L = laws.get("hybrid5")
    free = fit.fit_linear(L, pr, "mu", 2000, "none")
    con = fit.fit_linear(L, pr, "mu", 2000, "physical")
    assert not free.extra["physical"]
    assert con.extra["physical"] and con.active
    simp = fit.fit_linear(L, pr, "mu", 2000, "simplex")
    assert np.all(simp.w >= -1e-9) and simp.w.sum() <= 1 + 1e-9


def test_power2_recovers_exact():
    pr = Profile(MU, laws.power2(0.65, 0.6, MU))
    f = fit.fit_power2(pr, "r", 2000)
    assert np.allclose(f.w, [0.65, 0.6], atol=1e-8)


def test_quadrature_weights_integrate_measure():
    pr = law_profile("hybrid5", W5)
    for measure, total in (("mu", 1.0), ("r", 1.0), ("area", 1.0)):
        for rs in (None, 2000):
            _, _, q = quadrature(pr, measure, rs)
            assert np.sum(q) == pytest.approx(total, rel=1e-12)


def test_limb_cut_recovers_embedded_photosphere():
    # photosphere of radius 1 inside a shell of radius R; exponential tail beyond
    R, h = 1.004, 3e-4
    L = laws.get("hybrid4")
    w = np.array([0.1, 0.2, 0.15, 0.2])
    r_shell = np.sort(np.concatenate([np.linspace(0, 0.99, 50), np.linspace(0.99, 1.0, 80)[1:]]))
    r_true = r_shell * R
    mu_true = np.sqrt(np.clip(1 - r_true ** 2, 0, None))
    I = np.where(r_true <= 1, L.intensity(w, mu_true),
                 L.intensity(w, 0.0) * np.exp(-np.clip(r_true - 1, 0, None) / h))
    mu_shell = np.sqrt(1 - r_shell ** 2)
    pr, rc = limb_cut(mu_shell, I, "peak")
    assert rc == pytest.approx(1 / R, abs=2e-4)
    m = pr.mu > 0.05
    assert np.max(np.abs(pr.y[m] - L.intensity(w, pr.mu[m]) / L.intensity(w, 1.0))) < 5e-3


def test_light_curve_matches_squishierplanet_exact():
    L = laws.get("hybrid5")
    law = L.to_law(W5)
    for p, b in ((0.03, 0.2), (0.1, 0.6), (0.15, 0.9)):
        z = lc.chord_z(p, b, 40)
        F = lc.light_curve(lambda m: L.intensity(W5, m), z, p)
        Fex = np.array([sp.flux(law, X=zz, Y=0.0, a=p, b=p, psi=0.0) for zz in z])
        assert np.max(np.abs(F - Fex)) < 1e-12


def test_lc_errors_zero_for_identical_profiles():
    L = laws.get("hybrid4")
    f = lambda m: L.intensity(np.array([0.1, 0.2, 0.15, 0.2]), m)
    assert lc.lc_errors(f, f)["max_ppm"] == 0.0
