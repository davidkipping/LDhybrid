"""D14: physicality of coefficients summed from the coefficient spectra over
narrow top-hat bins (R = 100, 300, 1000; 0.6–5.3 µm), for ~150 PHOENIX models
spread over Teff, log g and [M/H]. The per-R=3000-bin fits are unconstrained,
so a sum over a narrow bin can leave the physical region.

    .venv/bin/python studies/rev/d14_narrow_bins.py     # writes d14.json
"""
import json

import numpy as np

import rc
from ldhybrid import fit, laws
from ldhybrid.paths import DATA
from ldhybrid.specint import SpecInt

SPEC = DATA / "products" / "spectra" / "phoenix"
RES = (100, 300, 1000)
LAWS = ("hybrid2", "hybrid4", "hybrid5")
X = np.linspace(0.0, 1.0, fit.N_CONSTRAINT)


def bins(R, lo=6000.0, hi=53000.0):
    n = int(np.floor(np.log(hi / lo) * R))
    return lo * np.exp(np.arange(n + 1) / R)


def main():
    files = sorted(SPEC.glob("*.npz"))
    rng = np.random.default_rng(3)
    pick = rng.choice(len(files), 150, replace=False)
    z0 = np.load(files[0])
    edges0 = z0["lam_edges"]
    dummy = SpecInt(edges0, np.array([0.0, 1.0]), np.zeros((len(edges0) - 1, 2)))
    Wm = {}
    for R in RES:
        e = bins(R)
        rows = []
        for a, b in zip(e[:-1], e[1:]):
            lam = np.array([a - 1e-3, a, b, b + 1e-3])
            rows.append(dummy.bin_throughput(lam, np.array([0.0, 1.0, 1.0, 0.0])))
        Wm[R] = (np.array(rows), e)
    G = {n: laws.get(n).dT(X) for n in LAWS}
    out = []
    for i in pick:
        z = np.load(files[i])
        assert np.array_equal(z["lam_edges"], edges0)
        teff, logg, mh = float(z["teff"]), float(z["logg"]), float(z["mh"])
        rec = dict(teff=teff, logg=logg, mh=mh, res={})
        for R, (M, e) in Wm.items():
            I1 = M @ z["I1"].astype(float)
            rec["res"][R] = {}
            for n in LAWS:
                w = (M @ z[f"a_{n}"].astype(float)) / I1[:, None]
                mono = np.all(-(w @ G[n].T) >= -1e-10, axis=1)
                pos = 1.0 - w.sum(1) >= -1e-10
                simplex = np.all(w >= -1e-10, axis=1) & pos
                bad = ~(mono & pos)
                rec["res"][R][n] = dict(n=int(len(w)), unphysical=int(bad.sum()),
                                        outside_simplex=int((~simplex).sum()),
                                        lam_bad=((e[:-1][bad] + e[1:][bad]) / 2).round(0).tolist())
        out.append(rec)
    (rc.HERE / "d14.json").write_text(json.dumps(out))
    print("done", len(out))


if __name__ == "__main__":
    main()
