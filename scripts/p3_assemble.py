"""P3: assemble chunk products into one HDF5 file per grid and long-format
gzipped CSV tables per grid × law × method.

    .venv/bin/python scripts/p3_assemble.py phoenix mps1 mps2 kurucz

Outputs (data/products/):
  <grid>.h5                       all arrays: keys (teff, logg, mh), bands,
                                  <law>|<method>|<field>
  tables/<grid>/<law>_<method>.csv.gz
      teff, logg, mh, band, r_c, w1..wn | (c, alpha) | (a1..a4) | (u1, u2),
      q1..qn and c0, c2, c4, d1..dK (hybrids), rms_I, max_abs_I, flux_err,
      lc_ppm, physical, bound
"""

import gzip
import sys

import h5py
import numpy as np

from ldhybrid import laws, qspace
from ldhybrid.paths import DATA
from ldhybrid.products import LAWS, METHODS

ROOT = DATA / "products"
PAR = {"hybrid2": ["w1", "w2"], "hybrid4": ["w1", "w2", "w3", "w4"],
       "hybrid5": ["w1", "w2", "w3", "w4", "w5"], "quadratic": ["u1", "u2"],
       "claret4": ["a1", "a2", "a3", "a4"], "power2": ["c", "alpha"]}


def load_chunks(grid):
    files = sorted((ROOT / "bands" / grid).glob("chunk_*.npz"))
    parts = [dict(np.load(f)) for f in files]
    out = {k: np.concatenate([p[k] for p in parts]) for k in parts[0] if k != "bands"}
    out["bands"] = parts[0]["bands"]
    return out


def write_h5(grid, A):
    with h5py.File(ROOT / f"{grid}.h5", "w") as h:
        keys = A["keys"]
        h["teff"] = keys[:, 0].astype(float)
        h["logg"] = keys[:, 1].astype(float)
        h["mh"] = keys[:, 2].astype(float)
        h["bands"] = A["bands"].astype("S")
        for k, v in A.items():
            if k not in ("keys", "bands"):
                h.create_dataset(k, data=v, compression="gzip")
        h.attrs["layout"] = "arrays indexed [model, band(, parameter)]"


def c_form(name, w):
    """Expanded coefficients: I/I(1) = c0 + c2 μ² + c4 μ⁴ + Σ d_k (μ² + ε_k)⁻²."""
    L = laws.get(name).to_law(w)
    c = np.zeros(3)
    c[: L.c.size] = L.c
    return list(c) + list(L.d)


def write_tables(grid, A):
    tdir = ROOT / "tables" / grid
    tdir.mkdir(parents=True, exist_ok=True)
    keys, bands = A["keys"], A["bands"]
    for name in LAWS:
        hyb = name.startswith("hybrid")
        for m in METHODS:
            k = f"{name}|{m}"
            W = A[k + "|w"]
            npar = W.shape[-1]
            head = ["teff", "logg", "mh", "band", "r_c"] + PAR[name]
            if hyb:
                K = npar - 2 if name != "hybrid2" else 1
                head += [f"q{j + 1}" for j in range(npar)] + ["c0", "c2", "c4"] + \
                        [f"d{j + 1}" for j in range(K)]
            head += ["rms_I", "max_abs_I", "flux_err", "lc_ppm", "physical", "bound"]
            path = tdir / f"{name}_{'r' if m == 'r' else 'mu'}.csv.gz"
            with gzip.open(path, "wt") as f:
                f.write(",".join(head) + "\n")
                for i in range(len(keys)):
                    t, g, mh = keys[i]
                    for j, b in enumerate(bands):
                        if not A["covered"][i, j]:
                            continue
                        w = W[i, j]
                        vals = [f"{float(t):.0f}", f"{float(g):.2f}", f"{float(mh):+.2f}", b,
                                f"{A['r_c'][i, j]:.6f}"] + [f"{x:.6f}" for x in w]
                        if hyb:
                            vals += [f"{x:.6f}" for x in qspace.to_q(name, w)]
                            vals += [f"{x:.8e}" for x in c_form(name, w)]
                        vals += [f"{A[k + '|rms_I'][i, j]:.3e}", f"{A[k + '|max_abs_I'][i, j]:.3e}",
                                 f"{A[k + '|flux_err'][i, j]:+.3e}", f"{A[k + '|lc_ppm'][i, j]:.3f}",
                                 str(int(A[k + "|physical"][i, j])), str(int(A[k + "|bound"][i, j]))]
                        f.write(",".join(vals) + "\n")
            print("  ", path.relative_to(ROOT), flush=True)


if __name__ == "__main__":
    for grid in sys.argv[1:]:
        A = load_chunks(grid)
        print(grid, len(A["keys"]), "models", len(A["bands"]), "bands", flush=True)
        write_h5(grid, A)
        write_tables(grid, A)
