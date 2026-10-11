"""Recompute the products of ONE law in the existing band chunks, in place
(e.g. after a change of that law's definition), without touching the others.

    .venv/bin/python scripts/p3_update_law.py hybrid4 phoenix mps1 mps2 [--workers N]

For every chunk in data/products/bands/<grid>/ the models are reloaded and the
law's arrays (both methods: w, rms_I, max_abs_I, flux_err, lc_ppm, physical,
bound) are recomputed exactly as products.model_products does and written back
atomically. Resumable: a chunk is marked done in <grid>/.<law>_updated.
"""

import argparse
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from p3_tables import OUT, load  # noqa: E402

from ldhybrid import bands as B  # noqa: E402
from ldhybrid import fit, laws, v1  # noqa: E402
from ldhybrid.products import METHODS, _cols, _grid, constraint_for  # noqa: E402
from ldhybrid.profile import band_profile  # noqa: E402


def typed_key(grid, k):
    t, g, mh = k
    return (int(float(t)), float(g), mh) if grid == "phoenix" else (int(float(t)), float(g), float(mh))


def update_chunk(grid, path, law):
    A = dict(np.load(path))
    bands = [str(b) for b in A["bands"]]
    g = _grid()
    for i, k in enumerate(A["keys"]):
        spec = load(grid, typed_key(grid, k))
        spherical = spec.meta.get("geometry") == "spherical"
        for j, b in enumerate(bands):
            if not A["covered"][i, j]:
                continue
            lam, S = B.response(b)
            pr = band_profile(spec, lam, S, limb=v1.LIMB) if spherical else band_profile(spec, lam, S)
            truth = g.columns(pr)
            for m in METHODS:
                key = f"{law}|{m}"
                f = fit.fit_linear(laws.get(law), pr, m, v1.RESAMPLE, constraint_for(law), v1.FCM)
                model = g.combine(*_cols(law), f.w)
                A[key + "|w"][i, j] = f.w
                A[key + "|rms_I"][i, j] = f.rms
                A[key + "|max_abs_I"][i, j] = f.max_abs
                A[key + "|flux_err"][i, j] = f.flux_err
                A[key + "|lc_ppm"][i, j] = g.errors(*model, *truth)["max_ppm"]
                A[key + "|physical"][i, j] = bool(f.extra["physical"])
                A[key + "|bound"][i, j] = bool(f.active)
    tmp = path.with_name(path.stem + ".tmp.npz")
    np.savez_compressed(tmp, **A)
    tmp.replace(path)
    return path.name


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("law")
    ap.add_argument("grids", nargs="+")
    ap.add_argument("--workers", type=int, default=6)
    a = ap.parse_args()
    for grid in a.grids:
        d = OUT / grid
        mark = d / f".{a.law}_updated"
        done = set(mark.read_text().split()) if mark.exists() else set()
        todo = [p for p in sorted(d.glob("chunk_*.npz")) if ".tmp" not in p.name and p.name not in done]
        print(f"{grid}: {len(todo)} chunks to update ({len(done)} done)", flush=True)
        t0 = time.time()
        with ProcessPoolExecutor(a.workers) as ex, open(mark, "a") as fm:
            futs = [ex.submit(update_chunk, grid, p, a.law) for p in todo]
            for n, f in enumerate(as_completed(futs), 1):
                fm.write(f.result() + "\n")
                fm.flush()
                if n % 20 == 0 or n == len(todo):
                    el = time.time() - t0
                    print(f"  {grid} {n}/{len(todo)}  {el / 60:.1f} min  eta {el / n * (len(todo) - n) / 60:.0f} min",
                          flush=True)
    print("ALL DONE", flush=True)


if __name__ == "__main__":
    sys.exit(main())
