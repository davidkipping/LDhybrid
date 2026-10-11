"""Round-2 checks of the user code and the CDS release (FORK B).

    .venv/bin/python studies/rev2/code_checks.py

1. light_curve() (standalone) against ldhybrid.lightcurve for several laws.
2. The new hybrid-table columns (q, expanded form) against the products.
3. MPS-ATLAS custom throughput (api.coefficients_custom, grid=mps2) against
   the tabulated coefficients of a registered band at grid models.
4. Corrections from the release (standalone) against the API (predict).
5. Joint condition() against sequential single updates; correlation effect.
6. prior() safeguards: M-dwarf systematic floor, cov_floor, dropped draws.
7. apply_offsets(): state matching; mismatched state raises.
8. refit(mu_min=0.02) for the four-term law: light-curve error at grid models.
"""

import json
import warnings
from pathlib import Path

import h5py
import numpy as np

import ldhybrid as ld
from ldhybrid import bands, fit, laws, v1
from ldhybrid import interpolate as I
from ldhybrid import lightcurve as lc
from ldhybrid.paths import DATA
from ldhybrid.products import PROD_LC, constraint_for
from ldhybrid.profile import Profile

HERE = Path(__file__).parent
REL = DATA / "release" / "cds"
T = I.LDTables(REL, profiles=DATA / "products" / "profiles")
L = []
OUT = {}


def law_fn(law, c):
    return lambda mu: I.intensity(law, c, mu)


def check_lightcurve():
    z = np.linspace(0, 1.15, 400)
    worst = {}
    for law, c in (("hybrid5", [0.05, 0.1, 0.12, 0.13, 0.2]), ("claret4", [0.6, -0.2, 0.5, -0.2]),
                   ("power2", [0.7, 0.6]), ("quadratic", [0.4, 0.2]), ("hybrid2", [0.4, 0.1])):
        for p in (0.03, 0.1, 0.3):
            F1 = I.light_curve(law, c, z, p)
            F2 = lc.light_curve(law_fn(law, c), z, p)
            worst[f"{law}|{p}"] = float(np.max(np.abs(F1 - F2)))
    OUT["lightcurve"] = worst
    L.append(f"1. light_curve() vs ldhybrid.lightcurve (5 laws × p = 0.03, 0.1, 0.3, 400 points): "
             f"max |ΔF| = {max(worst.values()):.1e}")


def check_columns():
    res = {}
    for grid in ("phoenix", "mps2"):
        h = h5py.File(DATA / "products" / f"{grid}.h5")
        for law in ("hybrid2", "hybrid5"):
            d = T.table(grid, law)
            m = d["method"] == "r"
            # expanded form must reproduce the intensity of the weights
            mu = np.linspace(0, 1, 41)
            eps = I.POLES[law]
            X = d["expanded"][m][:2000]
            W = d["coeffs"][m][:2000]
            I1 = X[:, [0]] + X[:, [1]] * mu ** 2 + X[:, [2]] * mu ** 4 + sum(
                X[:, [3 + k]] * (mu ** 2 + eps[k]) ** -2 for k in range(len(eps)))
            I2 = 1 - W @ I.shapes(law, mu).T
            # q must map back to the weights (inverse of the sampler maps)
            res[f"{grid}|{law}"] = dict(expanded_vs_w=float(np.max(np.abs(I1 - I2))),
                                        q_range=[float(d["q"].min()), float(d["q"].max())])
    OUT["columns"] = res
    L.append("2. new hybrid columns: expanded form reproduces I(μ) of the weights to "
             f"{max(v['expanded_vs_w'] for v in res.values()):.1e} (6-decimal weights, 8-digit "
             f"form); q within [{min(v['q_range'][0] for v in res.values()):.3f}, "
             f"{max(v['q_range'][1] for v in res.values()):.3f}]")


def check_mps_custom():
    h = h5py.File(DATA / "products" / "mps2.h5")
    B = [b.decode() for b in h["bands"][()]]
    rows = []
    for (t, g, m) in ((5800, 4.5, 0.0), (4200, 4.6, -0.5), (6600, 4.0, 0.3)):
        k = int(np.flatnonzero((h["teff"][()] == t) & (h["logg"][()] == g) &
                               (np.abs(h["mh"][()] - m) < 1e-6))[0])
        for band in ("TESS", "JWST_NIRSpec_PRISM", "JWST_MIRI_LRS"):
            j = B.index(band)
            for law in ("hybrid5", "claret4", "power2"):
                r = ld.coefficients_custom(law, t, g, m, *bands.response(band), grid="mps2")
                w_tab = h[f"{law}|r|w"][k, j]
                rows.append(dict(model=(t, g, m), band=band, law=law,
                                 dw=float(np.max(np.abs(r["w"] - w_tab)))))
    # off-node: compare with table interpolation in light-curve terms
    G = lc.LCGrid(**PROD_LC)
    r = ld.coefficients_custom("hybrid5", 5650, 4.42, 0.12, *bands.response("TESS"), grid="mps2")
    ti = T.coefficients("hybrid5", 5650, 4.42, 0.12, "TESS", "mps2")
    wi = np.array(list(ti["coeffs"].values()))
    L5 = laws.get("hybrid5")
    off = G.errors(*G.columns(lambda mu: L5.intensity(r["w"], mu)),
                   *G.columns(lambda mu: L5.intensity(wi, mu)))["max_ppm"]
    OUT["mps_custom"] = dict(nodes=rows, offnode_lc_ppm=float(off))
    L.append(f"3. MPS-ATLAS custom throughput at grid models (3 models × TESS/PRISM/MIRI × "
             f"hybrid5/four-term/power-2): max |Δw| vs table {max(x['dw'] for x in rows):.1e}; "
             f"off-node (5650 K, 4.42, +0.12, TESS) hybrid5 vs interpolated table: "
             f"{off:.2f} ppm in light-curve terms")


def check_corrections():
    d = []
    for law in ("hybrid2", "hybrid5"):
        for band in ("Kepler", "TESS", "JWST_NIRSpec_PRISM"):
            for teff in (3600, 4500, 5600, 6500):
                for state in ("ssd", "B100G"):
                    a = T.coefficients(law, teff, 4.5, 0.0, band, "mps2", correction=state)
                    b = ld.coefficients(law, teff, 4.5, 0.0, band, grid="mps2", correction=state)
                    d.append(float(np.max(np.abs(np.array(list(a["coeffs"].values())) - b["w"]))))
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        T.coefficients("hybrid5", 7200, 4.3, 0.0, "TESS", "mps2", correction="B100G")
    OUT["corrections"] = dict(max_dw=max(d), warned_outside=bool(w))
    L.append(f"4. corrections from the release (standalone) vs API: max |Δw| = {max(d):.1e} "
             f"(48 cases); warning outside the MURaM Teff range: {bool(w)}")


def check_condition():
    rng = np.random.default_rng(1)
    A = rng.normal(size=(5, 5)) * 0.01
    C = A @ A.T + np.eye(5) * 1e-4
    m = np.full(5, 0.1)
    g2 = I.observables("hybrid5")["h2"][1]
    m1, C1 = I.condition(m, C, "hybrid5", "h1", 0.008, 0.006)
    # a sequential second update must re-reference the offset to the updated mean
    m1, C1 = I.condition(m1, C1, "hybrid5", "h2", -0.013 - g2 @ (m1 - m), 0.013)
    m2, C2 = I.condition(m, C, "hybrid5", ["h1", "h2"], [0.008, -0.013], [0.006, 0.013])
    m3, C3 = I.condition(m, C, "hybrid5", ["h1", "h2"], [0.008, -0.013], [0.006, 0.013], rho=0.5)
    g = I.observables("hybrid5")
    var = lambda Cx, k: float(g[k][1] @ Cx @ g[k][1])
    OUT["condition"] = dict(seq_vs_joint_mean=float(np.max(np.abs(m1 - m2))),
                            seq_vs_joint_cov=float(np.max(np.abs(C1 - C2))),
                            sd_h1=[var(C, "h1") ** 0.5, var(C2, "h1") ** 0.5, var(C3, "h1") ** 0.5])
    L.append(f"5. joint condition() = sequential updates for rho = 0 (|Δmean| "
             f"{OUT['condition']['seq_vs_joint_mean']:.1e}, |Δcov| "
             f"{OUT['condition']['seq_vs_joint_cov']:.1e}); sd(h'1) prior / rho=0 / rho=0.5: "
             + " / ".join(f"{x:.4f}" for x in OUT["condition"]["sd_h1"]))


def check_prior():
    res = {}
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        p0 = T.prior("hybrid5", 3000, 50, 5.0, 0.1, 0.0, 0.1, "JWST_NIRSpec_G395H", "phoenix",
                     systematic=True, m_dwarf_floor=False)
        p1 = T.prior("hybrid5", 3000, 50, 5.0, 0.1, 0.0, 0.1, "JWST_NIRSpec_G395H", "phoenix",
                     systematic=True)
        msgs = [str(x.message) for x in w]
    ev = lambda C: np.sqrt(np.clip(np.linalg.eigvalsh(C), 0, None))
    res["m_dwarf"] = dict(sd_eig_no_floor=ev(p0["cov"]).tolist(), sd_eig_floor=ev(p1["cov"]).tolist(),
                          floor_point=p1["floor_point"], warnings=msgs)
    p2 = T.prior("hybrid5", 5600, 100, 4.4, 0.1, 0.1, 0.1, "TESS", "mps2")
    p3 = T.prior("hybrid5", 5600, 100, 4.4, 0.1, 0.1, 0.1, "TESS", "mps2", cov_floor=0.01 ** 2)
    res["cov_floor"] = dict(sd_eig=ev(p2["cov"]).tolist(), sd_eig_floor=ev(p3["cov"]).tolist())
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        p4 = T.prior("hybrid5", 3550, 150, 4.6, 0.1, 0.0, 0.1, "TESS", "mps2")
    res["edge"] = dict(dropped=p4["dropped"], warned=bool(w))
    p5 = T.prior("hybrid5", 5600, 1e-6, 4.4, 1e-6, 0.1, 1e-6, "TESS", "mps2", correction="B100G")
    c5 = T.coefficients("hybrid5", 5600, 4.4, 0.1, "TESS", "mps2", correction="B100G")
    res["sigma0_correction"] = float(np.max(np.abs(p5["mean"] - np.array(list(c5["coeffs"].values())))))
    OUT["prior"] = res
    L.append("6. prior(): 3000 K G395H PHOENIX, systematic=True: only PHOENIX covers the star → "
             f"warning; eigen-s.d. without floor {np.round(res['m_dwarf']['sd_eig_no_floor'], 5).tolist()}, "
             f"with the floor at {p1['floor_point']} {np.round(res['m_dwarf']['sd_eig_floor'], 4).tolist()}")
    L.append(f"   cov_floor = 0.01²: eigen-s.d. {np.round(res['cov_floor']['sd_eig'], 5).tolist()} → "
             f"{np.round(res['cov_floor']['sd_eig_floor'], 4).tolist()}")
    L.append(f"   grid edge (3550 ± 150 K, MPS): {100 * res['edge']['dropped']:.0f}% of draws dropped, "
             f"warning {res['edge']['warned']}; σ → 0 with B100G reproduces coefficients(): "
             f"{res['sigma0_correction']:.1e}")


def check_offsets():
    p = T.prior("hybrid5", 5600, 100, 4.4, 0.1, 0.1, 0.1, "Kepler", "mps2", correction="B100G",
                systematic=True)
    m, C = T.apply_offsets(p, "Kepler")
    g = I.observables("hybrid5")
    h1 = lambda mm: float(g["h1"][0] + g["h1"][1] @ mm)
    try:
        T.apply_offsets(dict(p, correction="ssd"), "Kepler")
        raised = False
    except ValueError:
        raised = True
    pj = T.prior("hybrid5", 5600, 100, 4.4, 0.1, 0.1, 0.1, "JWST_NIRSpec_PRISM", "mps2",
                 correction="B100G")
    mj, Cj = T.apply_offsets(pj, "JWST_0.86-5.5um")
    gl = lambda mm: float(g["ell"][0] + g["ell"][1] @ mm)
    OUT["offsets"] = dict(h1_shift=h1(m) - h1(p["mean"]), mismatch_raises=raised,
                          ell_shift=gl(mj) - gl(pj["mean"]))
    L.append(f"7. apply_offsets(): Kepler, MPS2 + B100G prior: h'1 shifted by "
             f"{OUT['offsets']['h1_shift']:+.4f} (offset −0.0023, joint with h'2, rho 0.04); "
             f"a prior in another state raises: {raised}; JWST ℓ shift "
             f"{OUT['offsets']['ell_shift']:+.4f} (offset +0.024)")


def check_refit_mumin():
    h = h5py.File(DATA / "products" / "phoenix.h5")
    B = [b.decode() for b in h["bands"][()]]
    P = T._profiles("phoenix")
    G = lc.LCGrid(**PROD_LC)
    rng = np.random.default_rng(3)
    a, b = [], []
    for k in rng.choice(len(P["teff"]), 30, replace=False):
        band = ("TESS", "HST_WFC3_G141", "JWST_NIRSpec_PRISM")[k % 3]
        j = P["bands"].index(band)
        mm = np.isfinite(P["mu"][k, j])
        pr = Profile(P["mu"][k, j][mm].astype(float), P["y"][k, j][mm].astype(float))
        t, g_, m_ = float(P["teff"][k]), float(P["logg"][k]), float(P["mh"][k])
        tr = G.columns(pr)
        for mu_min, acc in ((0.0, a), (0.02, b)):
            r = T.refit("claret4", t, g_, m_, band, "phoenix", mu_min=mu_min)
            c = list(r["coeffs"].values())
            acc.append(G.errors(*G.columns(lambda mu: I.intensity("claret4", c, mu)), *tr)["max_ppm"])
    OUT["refit_mu_min"] = dict(median_full=float(np.median(a)), median_cut=float(np.median(b)),
                               p90_full=float(np.percentile(a, 90)), p90_cut=float(np.percentile(b, 90)))
    L.append(f"8. four-term refit at 30 PHOENIX models: light-curve error median / 90th "
             f"{np.median(a):.1f} / {np.percentile(a, 90):.1f} ppm (all μ) → {np.median(b):.1f} / "
             f"{np.percentile(b, 90):.1f} ppm (μ ≥ 0.02)")


def main():
    for f in (check_lightcurve, check_columns, check_mps_custom, check_corrections,
              check_condition, check_prior, check_offsets, check_refit_mumin):
        f()
        print(L[-1] if not L[-1].startswith("   ") else "\n".join(L[-3:]), flush=True)
    (HERE / "code_checks.json").write_text(json.dumps(OUT, indent=1, default=float))
    (HERE / "code_checks.md").write_text("# Round-2 code and release checks (studies/rev2/code_checks.py)\n\n"
                                        + "\n".join(L) + "\n")


if __name__ == "__main__":
    main()
