"""C6: how the linear-interpolation error scales with the grid spacing.

PHOENIX (Teff step h = 100 K): at a node T, interpolate the tabulated
coefficients linearly across a gap of 2h (T ± 100) and of 3h (T − 100 and
T + 200, weights 2/3 and 1/3; and T − 200 and T + 100), and score each
against the light curve of the tabulated (post-cut) profile at T. Fitting
error ∝ gap^k gives k and the extrapolated error at the native spacing.
MPS-ATLAS Set 2 in [M/H] (native step 0.05 dex between −1 and +0.5): gaps of
0.1, 0.2 and 0.3 dex.
Laws: hybrid5, four-term, power-2 (r-weighted). Production metric.

    .venv/bin/python studies/rev2/c6_interp.py
"""
import json
from pathlib import Path

import h5py
import numpy as np

import rc
from ldhybrid import laws
from ldhybrid.profile import Profile

ROOT = Path(__file__).resolve().parents[2]
P = ROOT / "data" / "products"
LAWS = ["hybrid5", "claret4", "power2"]


def law_fn(name, w):
    if name == "power2":
        return lambda mu: laws.power2(w[0], w[1], mu)
    L = laws.get(name)
    return lambda mu: L.intensity(w, mu)


_H = {}


def data(grid):
    if grid not in _H:
        H = h5py.File(P / f"{grid}.h5")
        Pr = h5py.File(P / "profiles" / f"{grid}.h5")
        keys = {(float(t), float(g), round(float(m), 3)): i
                for i, (t, g, m) in enumerate(zip(H["teff"][()], H["logg"][()], H["mh"][()]))}
        pk = {(float(t), float(g), round(float(m), 3)): i
              for i, (t, g, m) in enumerate(zip(Pr["teff"][()], Pr["logg"][()], Pr["mh"][()]))}
        _H[grid] = dict(H=H, P=Pr, keys=keys, pk=pk, hb=[b.decode() for b in H["bands"][()]],
                        pb=[b.decode() for b in Pr["bands"][()]],
                        W={n: H[f"{n}|r|w"][()] for n in LAWS}, cov=H["covered"][()])
    return _H[grid]


def task(job):
    grid, k0, band, schemes = job
    D = data(grid)
    jh, jp = D["hb"].index(band), D["pb"].index(band)
    ip = D["pk"][k0]
    m = np.isfinite(D["P"]["mu"][ip, jp])
    truth = Profile(D["P"]["mu"][ip, jp][m].astype(float), D["P"]["y"][ip, jp][m].astype(float))
    g = rc.lcgrid()
    tr = g.columns(truth)
    out = dict(grid=grid, key=k0, band=band)
    tab = {n: g.columns(law_fn(n, D["W"][n][D["keys"][k0], jh])) for n in LAWS}
    for n in LAWS:
        out[f"tab|{n}"] = g.errors(*tab[n], *tr)["max_ppm"]
    for sname, corners in schemes.items():
        for n in LAWS:
            w = sum(wt * D["W"][n][D["keys"][c], jh] for c, wt in corners)
            cm = g.columns(law_fn(n, w))
            out[f"{sname}|{n}"] = g.errors(*cm, *tr)["max_ppm"]
            out[f"{sname}|{n}|vs_tab"] = g.errors(*cm, *tab[n])["max_ppm"]
    return out


def jobs_phoenix(n=1500, seed=0):
    D = data("phoenix")
    rng = np.random.default_rng(seed)
    keys = list(D["keys"])
    out = []
    while len(out) < n:
        t, g, mh = keys[rng.integers(len(keys))]
        need = [(t + d, g, mh) for d in (-200, -100, 100, 200)]
        if not all(x in D["keys"] and x in D["pk"] for x in need) or (t, g, mh) not in D["pk"]:
            continue
        band = D["hb"][rng.integers(len(D["hb"]))]
        j = D["hb"].index(band)
        if not all(D["cov"][D["keys"][x], j] for x in need + [(t, g, mh)]):
            continue
        sch = {"2h": [((t - 100, g, mh), 0.5), ((t + 100, g, mh), 0.5)],
               "3h_a": [((t - 100, g, mh), 2 / 3), ((t + 200, g, mh), 1 / 3)],
               "3h_b": [((t - 200, g, mh), 1 / 3), ((t + 100, g, mh), 2 / 3)]}
        out.append(("phoenix", (t, g, mh), band, sch))
    return out


def jobs_mps(n=600, seed=1):
    D = data("mps2")
    rng = np.random.default_rng(seed)
    keys = [k for k in D["keys"] if -0.9 <= k[2] <= 0.4]
    out = []
    while len(out) < n:
        t, g, mh = keys[rng.integers(len(keys))]
        r3 = lambda x: round(x, 3)
        need = [(t, g, r3(mh + d)) for d in (-0.1, -0.05, 0.05, 0.1, 0.2)]
        if not all(x in D["keys"] and x in D["pk"] for x in need) or (t, g, mh) not in D["pk"]:
            continue
        band = D["hb"][rng.integers(len(D["hb"]))]
        sch = {"0.1": [((t, g, r3(mh - 0.05)), 0.5), ((t, g, r3(mh + 0.05)), 0.5)],
               "0.2": [((t, g, r3(mh - 0.1)), 0.5), ((t, g, r3(mh + 0.1)), 0.5)],
               "0.3": [((t, g, r3(mh - 0.1)), 2 / 3), ((t, g, r3(mh + 0.2)), 1 / 3)]}
        out.append(("mps2", (t, g, mh), band, sch))
    return out


def main():
    jobs = jobs_phoenix() + jobs_mps()
    rows = rc.run_parallel(task, jobs)
    (rc.HERE / "c6_interp.json").write_text(json.dumps(rows, default=str))
    L = ["## C6 — interpolation error vs grid spacing (production metric, ppm)", ""]
    ph = [r for r in rows if r["grid"] == "phoenix"]
    L += [f"PHOENIX Teff (native step h = 100 K), {len(ph)} node-bands:", "",
          "| law | gap 2h (200 K) med / 90th | gap 3h (300 K) med / 90th | k (median ratio) | native h (extrapolated) med / 90th |",
          "|---|---|---|---|---|"]
    for n in LAWS:
        e2 = np.array([r[f"2h|{n}"] for r in ph])
        e3 = np.array([0.5 * (r[f"3h_a|{n}"] + r[f"3h_b|{n}"]) for r in ph])
        k = np.log(np.median(e3) / np.median(e2)) / np.log(1.5)
        k90 = np.log(np.percentile(e3, 90) / np.percentile(e2, 90)) / np.log(1.5)
        L.append(f"| {n} | {np.median(e2):.1f} / {np.percentile(e2, 90):.1f} | {np.median(e3):.1f} / "
                 f"{np.percentile(e3, 90):.1f} | {k:.2f} (90th: {k90:.2f}) | "
                 f"{np.median(e2) * 0.5 ** k:.1f} / {np.percentile(e2, 90) * 0.5 ** k90:.1f} |")
    ms = [r for r in rows if r["grid"] == "mps2"]
    L += ["", f"MPS-ATLAS Set 2 [M/H] (native step 0.05 dex in −1…+0.5), {len(ms)} node-bands:", "",
          "| law | gap 0.1 dex | gap 0.2 dex | gap 0.3 dex | k (0.1→0.2) |", "|---|---|---|---|---|"]
    for n in LAWS:
        e = {s: np.array([r[f"{s}|{n}"] for r in ms]) for s in ("0.1", "0.2", "0.3")}
        k = np.log(np.median(e["0.2"]) / np.median(e["0.1"])) / np.log(2)
        L.append(f"| {n} | " + " | ".join(f"{np.median(e[s]):.2f} / {np.percentile(e[s], 90):.1f}" for s in e)
                 + f" | {k:.2f} |")
    L += ["", "Pure interpolation error (interpolated vs tabulated coefficients at the node), median / 90th:", "",
          "| law | grid | " + " | ".join(["gap 2h/0.1", "gap 3h/0.2", "(MPS 0.3)", "k", "native (extrap.)"]) + " |",
          "|---|---|---|---|---|---|---|"]
    for n in LAWS:
        e2 = np.array([r[f"2h|{n}|vs_tab"] for r in ph])
        e3 = np.array([0.5 * (r[f"3h_a|{n}|vs_tab"] + r[f"3h_b|{n}|vs_tab"]) for r in ph])
        k = np.log(np.median(e3) / np.median(e2)) / np.log(1.5)
        k90 = np.log(np.percentile(e3, 90) / np.percentile(e2, 90)) / np.log(1.5)
        L.append(f"| {n} | PHOENIX Teff | {np.median(e2):.1f} / {np.percentile(e2, 90):.1f} | "
                 f"{np.median(e3):.1f} / {np.percentile(e3, 90):.1f} | – | {k:.2f} | "
                 f"{np.median(e2) * 0.5 ** k:.1f} / {np.percentile(e2, 90) * 0.5 ** k90:.1f} |")
    for n in LAWS:
        e = {s_: np.array([r[f"{s_}|{n}|vs_tab"] for r in ms]) for s_ in ("0.1", "0.2", "0.3")}
        k = np.log(np.median(e["0.2"]) / np.median(e["0.1"])) / np.log(2)
        L.append(f"| {n} | MPS [M/H] | " + " | ".join(f"{np.median(e[s_]):.2f} / {np.percentile(e[s_], 90):.1f}"
                                                     for s_ in e) + f" | {k:.2f} | (native = 0.1 gap) |")
    tabs = {n: np.median([r[f"tab|{n}"] for r in ph]) for n in LAWS}
    L += ["", "Tabulated law error at the same PHOENIX nodes (median): " + ", ".join(f"{n} {v:.1f}" for n, v in tabs.items())]
    (rc.HERE / "c6_interp.md").write_text("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
