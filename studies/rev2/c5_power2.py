"""C5: power-2 physicality by Teff and band (rebuilt products), the c
distribution of unphysical rows, the T3 rows, and bound fractions of the
uniform-μ fits.

    .venv/bin/python studies/rev2/c5_power2.py      # writes c5_power2.json/.md
"""
import json
from pathlib import Path

import h5py
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).parent
P = ROOT / "data" / "products"
LAWS = ["hybrid5", "hybrid4", "hybrid2", "claret4", "quadratic"]
TBINS = [(2300, 3000), (3000, 3500), (3500, 4000), (4000, 99999)]


def main():
    out, L = {}, ["## C5 — power-2 physicality and uniform-μ constraint fractions", ""]
    for grid in ("phoenix", "mps1", "mps2"):
        h = h5py.File(P / f"{grid}.h5")
        cov = h["covered"][()]
        T = np.broadcast_to(h["teff"][()][:, None], cov.shape)
        B = [b.decode() for b in h["bands"][()]]
        g = out[grid] = {}
        for m in ("r", "mu"):
            W = h[f"power2|{m}|w"][()]
            c = W[..., 0]
            unph = cov & ~h[f"power2|{m}|physical"][()]
            g[f"unphys_{m}"] = float(unph.sum() / cov.sum())
            g[f"unphys_by_teff_{m}"] = {f"{a}-{b}": float(unph[(T >= a) & (T < b)].sum() /
                                                       max(1, cov[(T >= a) & (T < b)].sum()))
                                        for a, b in TBINS}
            cu = c[unph]
            if cu.size:
                g[f"c_unphys_{m}"] = dict(median=float(np.median(cu)), p90=float(np.percentile(cu, 90)),
                                          max=float(cu.max()), frac_gt2=float(np.mean(cu > 2)),
                                          n=int(cu.size))
            cool = T < 3500
            if cool.any():
                g[f"unphys_lt3500_by_band_{m}"] = {
                    b: float(unph[:, j][cool[:, j]].sum() / max(1, cov[:, j][cool[:, j]].sum()))
                    for j, b in enumerate(B)}
            g[f"bound_{m}"] = {law: float((cov & h[f"{law}|{m}|bound"][()]).sum() / cov.sum())
                               for law in LAWS}
    # T3 rows
    d = json.loads((ROOT / "studies" / "lc" / "lc_section.json").read_text())
    t3 = {k: v for k, v in d["coeffs"].items() if k.endswith("|power2")}
    out["t3_power2"] = {k.rsplit("|", 1)[0]: dict(c=v[0], alpha=v[1], physical=bool(0 <= v[0] <= 1))
                        for k, v in t3.items()}
    (HERE / "c5_power2.json").write_text(json.dumps(out, indent=1))

    for grid in ("phoenix", "mps1", "mps2"):
        g = out[grid]
        L.append(f"**{grid}**: power-2 unphysical r {100 * g['unphys_r']:.1f}%, μ {100 * g['unphys_mu']:.1f}%")
        for m in ("r", "mu"):
            L.append(f"* by Teff ({m}): " + ", ".join(f"{k} K {100 * v:.0f}%" for k, v in g[f'unphys_by_teff_{m}'].items()))
            if f"c_unphys_{m}" in g:
                cu = g[f"c_unphys_{m}"]
                L.append(f"* c of unphysical rows ({m}): median {cu['median']:.3f}, 90th {cu['p90']:.2f}, "
                         f"max {cu['max']:.1f}, c > 2 in {100 * cu['frac_gt2']:.1f}% (n = {cu['n']})")
            if f"unphys_lt3500_by_band_{m}" in g and m == "r":
                bb = g[f"unphys_lt3500_by_band_{m}"]
                top = sorted(bb.items(), key=lambda kv: -kv[1])
                L.append("* < 3500 K by band (r): " + ", ".join(f"{b} {100 * v:.0f}%" for b, v in top))
            L.append(f"* bound fractions ({m}): " + ", ".join(f"{k} {100 * v:.2f}%" for k, v in g[f'bound_{m}'].items()))
        L.append("")
    L.append("T3 power-2 (tabulated, r): " + "; ".join(
        f"{k} c={v['c']:.3f}{'' if v['physical'] else ' UNPHYSICAL'}" for k, v in out["t3_power2"].items()))
    (HERE / "c5_power2.md").write_text("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
