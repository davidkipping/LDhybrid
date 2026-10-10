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
