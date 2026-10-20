"""The standalone table interpolator (ldhybrid/interpolate.py)."""

import gzip

import numpy as np
import pytest

from ldhybrid import laws
from ldhybrid.interpolate import COEFFS, LDTables, intensity, is_physical, shapes
from ldhybrid.paths import DATA


def test_shapes_match_package_laws():
    mu = np.linspace(0, 1, 201)
    for law in ("hybrid2", "hybrid4", "hybrid5", "quadratic", "claret4"):
        w = np.random.default_rng(0).dirichlet(np.ones(len(COEFFS[law]) + 1))[:-1]
        assert np.allclose(intensity(law, w, mu), laws.get(law).intensity(w, mu), atol=1e-12)
    assert np.allclose(shapes("hybrid5", np.array([0.0, 1.0])),
                       [[1, 1, 1, 1, 1], [0, 0, 0, 0, 0]], atol=1e-12)


def _write_table(path, law):
    """A tiny synthetic table in the CDS layout: coefficients linear in the parameters."""
    rows = []
    for t in (5000, 5100, 5200):
        for g in (4.0, 4.5):
            for m in (-0.5, 0.0):
                c = 0.1 + 1e-4 * (t - 5000) + 0.02 * g + 0.01 * m + 0.01 * np.arange(len(COEFFS[law]))
                rows.append(f"{t:5d} {g:4.2f} {m:5.2f} {'TESS':25s} r  0.050000 " +
                            " ".join(f"{x:10.6f}" for x in c) +
                            " 1.000E-04 1.000E-03  1.000E-05     1.00 1 0")
    with gzip.open(path, "wt") as f:
        f.write("\n".join(rows) + "\n")


def test_exact_for_linear_coefficients(tmp_path):
    _write_table(tmp_path / "phoenix_hybrid5.dat.gz", "hybrid5")
    t = LDTables(tmp_path)
    r = t.coefficients("hybrid5", 5137, 4.31, -0.21, "TESS", grid="phoenix")
    exp = 0.1 + 1e-4 * 137 + 0.02 * 4.31 + 0.01 * -0.21 + 0.01 * np.arange(5)
    assert np.allclose(list(r["coeffs"].values()), exp, atol=1e-9) and r["inside"]
    with pytest.raises(ValueError):
        t.coefficients("hybrid5", 5400, 4.3, 0.0, "TESS", grid="phoenix")
    assert (tmp_path / "phoenix_hybrid5.npz").exists()               # cache written
    r2 = LDTables(tmp_path).coefficients("hybrid5", 5137, 4.31, -0.21, "TESS", grid="phoenix")
    assert r2["coeffs"] == r["coeffs"]


def test_is_physical():
    assert is_physical("hybrid5", [0.1, 0.2, 0.1, 0.1, 0.2])
    assert not is_physical("hybrid5", [0.5, 0.5, 0.3, 0.0, 0.0])     # limb intensity < 0
    assert not is_physical("power2", [1.2, 0.5])


CDS = DATA / "release" / "cds"


@pytest.mark.data
@pytest.mark.skipif(not (CDS / "phoenix_hybrid4.dat.gz").exists(), reason="CDS tables not built")
def test_agrees_with_package_api():
    import ldhybrid as ld
    t = LDTables(CDS, cache=False)
    a = t.coefficients("hybrid4", 4321, 4.37, -0.13, "HST_WFC3_G141", grid="phoenix")
    b = ld.coefficients("hybrid4", 4321, 4.37, -0.13, "HST_WFC3_G141", grid="phoenix")
    assert np.allclose(list(a["coeffs"].values()), b["w"], atol=2e-6)   # 6-decimal tables


def test_mu_cri_key_and_eigenvalue_floor(tmp_path):
    _write_table(tmp_path / "phoenix_hybrid5.dat.gz", "hybrid5")
    t = LDTables(tmp_path, cache=False)
    assert t.coefficients("hybrid5", 5100, 4.2, -0.2, "TESS", grid="phoenix")["mu_cri"] == 0.05
    # a Teff-only spread is rank one; the floor must lift every other direction
    p0 = t.prior("hybrid5", 5100, 30, 4.2, 0.0, -0.2, 0.0, "TESS", grid="phoenix")
    p = t.prior("hybrid5", 5100, 30, 4.2, 0.0, -0.2, 0.0, "TESS", grid="phoenix", cov_floor=1e-6)
    lam0, lam = np.linalg.eigvalsh(p0["cov"]), np.linalg.eigvalsh(p["cov"])
    assert lam0[0] < 1e-12 and np.isclose(lam[0], 1e-6) and lam0[-1] > 1e-5
    assert np.isclose(lam[-1], lam0[-1])                          # wide direction untouched


@pytest.mark.data
@pytest.mark.skipif(not (CDS / "offsets.dat").exists(), reason="CDS tables not built")
def test_offsets_check_the_passband():
    t = LDTables(CDS, cache=False)
    p = t.prior("hybrid5", 5600, 50, 4.4, 0.1, 0.0, 0.1, "TESS", grid="mps2", n=200)
    t.apply_offsets(p, "TESS")
    with pytest.raises(ValueError, match="not measured in TESS"):
        t.apply_offsets(p, "Kepler")
    t.apply_offsets(p, "Kepler", force=True)


@pytest.mark.data
@pytest.mark.skipif(not (CDS / "mps2_hybrid5.dat.gz").exists(), reason="CDS tables not built")
@pytest.mark.parametrize("law", ["hybrid4", "hybrid5"])
def test_standalone_files_reproduce_the_tables(law):
    import importlib.util
    from pathlib import Path
    spec = importlib.util.spec_from_file_location(
        "sa", Path(__file__).resolve().parents[1] / "standalone" / f"ldhybrid_{law}_mps2_tess.py")
    sa = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(sa)
    t = LDTables(CDS, cache=False)
    for x in [(5600, 4.4, 0.1), (4321, 4.37, -0.13), (6789, 4.1, -0.6)]:
        a = t.coefficients(law, *x, "TESS", grid="mps2", correction="B100G")
        b = sa.coefficients(*x, correction="B100G")
        assert np.allclose(list(a["coeffs"].values()), b["w"], atol=1e-7)   # float32 storage
        assert np.allclose(a["q"], b["q"], atol=1e-7) and a["physical"] == b["physical"]
    q = np.random.default_rng(0).uniform(size=(100, len(COEFFS[law])))
    w = sa.from_q(q)
    assert np.allclose([sa.to_q(v) for v in w], q, atol=1e-12)
    assert not any(np.all((np.array(sa.to_q(-v)) >= 0) & (np.array(sa.to_q(-v)) <= 1)) for v in w)


def test_grid_simplex_exact_continuous_and_filled():
    from ldhybrid.interpolate import GridSimplex
    rng = np.random.default_rng(0)
    axes = [np.array([0.0, 1.0, 3.0, 4.0]), np.array([0.0, 0.5, 1.5]), np.array([-1.0, 0.0, 2.0])]
    X = np.array(np.meshgrid(*axes, indexing="ij")).reshape(3, -1).T
    keep = ~np.all(X == [1.0, 0.5, 0.0], axis=1)        # one interior hole, filled along an axis
    keep &= ~np.all(X == [4.0, 1.5, 2.0], axis=1)       # a corner hole: no bracket on any axis
    A = rng.normal(size=(3, 2))
    G = GridSimplex(X[keep], np.column_stack([X[keep] @ A, rng.normal(size=keep.sum())]))
    assert (G.n_filled, G.n_empty) == (1, 1)
    assert np.allclose(G(X[keep]), G.Y)                               # exact at the models
    P = rng.uniform([0, 0, -1], [3.0, 1.4, 1.9], size=(2000, 3))     # away from the empty corner
    assert np.allclose(G(P)[:, :2], P @ A, atol=1e-12)                # exact for linear data
    D = rng.normal(size=P.shape)
    assert np.abs(G(P + 1e-9 * D) - G(P - 1e-9 * D)).max() < 1e-6    # continuous (no jumps)
    assert np.isnan(G([[3.9, 1.4, 1.9]])).all()                       # uses the empty corner
    assert np.isnan(G([[4.5, 0.2, 0.0]])).all()                       # outside the axes
    mw = G.model_weights([1.0, 0.5, 0.0])                             # the filled node: two models
    assert len(mw) == 2 and np.isclose(sum(mw.values()), 1.0)


@pytest.mark.data
@pytest.mark.skipif(not (CDS / "mps2_hybrid4.dat.gz").exists(), reason="CDS tables not built")
def test_tables_interpolant_continuous_at_grid_midpoints():
    # a Delaunay triangulation of these grids jumped by up to 1e-2 across such faces
    t = LDTables(CDS, cache=False)
    rng = np.random.default_rng(1)
    for grid in ("phoenix", "mps2"):
        f, _ = t._interpolator(grid, "hybrid4", "TESS", "r", None)
        X = f.X
        P = 0.5 * (X[rng.integers(len(X), size=2000)] + X[rng.integers(len(X), size=2000)])
        D = rng.normal(size=P.shape)
        a, b = f(P + 1e-7 * D), f(P - 1e-7 * D)
        ok = np.isfinite(a[:, 0]) & np.isfinite(b[:, 0])
        assert ok.mean() > 0.5 and np.abs(a - b)[ok, :4].max() < 1e-5


@pytest.mark.data
@pytest.mark.skipif(not (CDS / "mps2_hybrid5.dat.gz").exists(), reason="CDS tables not built")
@pytest.mark.parametrize("law", ["hybrid4", "hybrid5"])
def test_mlx_standalone_matches_numpy_and_differentiates(law):
    mx = pytest.importorskip("mlx.core")
    import importlib.util
    from pathlib import Path

    def load(name):
        spec = importlib.util.spec_from_file_location(
            name, Path(__file__).resolve().parents[1] / "standalone" / f"{name}.py")
        m = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(m)
        return m
    sa, sm = load(f"ldhybrid_{law}_mps2_tess"), load(f"ldhybrid_{law}_mps2_tess_mlx")
    rng = np.random.default_rng(2)
    P = rng.uniform([3400, 3.4, -2.1], [7600, 5.1, 1.1], size=(300, 3))
    r = sm.coefficients(P[:, 0], P[:, 1], P[:, 2], correction="B100G")   # float64
    w, inside = np.array(r["w"]), np.array(r["inside"])
    assert not inside.all() and np.isnan(w[~inside]).all()
    ref = np.array([sa.coefficients(*p, correction="B100G")["w"] for p in P[inside]])
    assert np.abs(w[inside] - ref).max() < 1e-14
    assert np.allclose(np.array(r["q"])[inside], [sa.to_q(v) for v in ref], atol=1e-12)
    # float32 gradient in Teff against float64 central differences
    T, g, m = rng.uniform(4000, 7000, 64), 4.37, -0.13
    c = rng.normal(size=sm.NW)
    gT = mx.grad(lambda t: (sm.coefficients(t, g, m)["w"] * mx.array(c.astype(np.float32))).sum())(
        mx.array(T.astype(np.float32)))
    fd = (np.array(sm.coefficients(T + 1e-3, g, m)["w"]) - np.array(sm.coefficients(T - 1e-3, g, m)["w"])) @ c / 2e-3
    assert np.median(np.abs(np.array(gT) - fd) / np.abs(fd)) < 1e-3
