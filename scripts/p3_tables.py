"""P3: compute band-table products for a grid, in resumable chunks.

    .venv/bin/python scripts/p3_tables.py phoenix|mps1|mps2|kurucz [--workers N]

Writes data/products/bands/<grid>/chunk_NNNN.npz (arrays over models × bands);
assemble with scripts/p3_assemble.py.
"""

import argparse
import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np

from ldhybrid import bands
from ldhybrid.paths import DATA, REDUCED

OUT = DATA / "products" / "bands"
CHUNK = 40
#: MIRI LRS only from Kurucz (v1.0 scope: FGK, 3500–6500 K)
BANDS = {"phoenix": [b for b in bands.BANDS if b != "JWST_MIRI_LRS"],
         "mps1": [b for b in bands.BANDS if b != "JWST_MIRI_LRS"],
         "mps2": [b for b in bands.BANDS if b != "JWST_MIRI_LRS"],
         "kurucz": ["JWST_MIRI_LRS"]}
SCOPE = dict(teff=(2300, 7500), logg=(3.5, 5.5), mh=(-2.0, 1.0))


def in_scope(teff, logg, mh):
    return (SCOPE["teff"][0] <= teff <= SCOPE["teff"][1]
            and SCOPE["logg"][0] <= logg <= SCOPE["logg"][1]
            and SCOPE["mh"][0] - 1e-9 <= mh <= SCOPE["mh"][1] + 1e-9)


def models(grid):
    if grid == "phoenix":
        cat = json.load(open(REDUCED / "phoenix" / "catalogue.json"))
        return sorted((e["teff"], e["logg"], e["mhs"]) for e in cat)
    if grid in ("mps1", "mps2"):
        from ldhybrid.grids import mpsatlas
        return [(t, g, mh) for mh, t, g in mpsatlas.catalogue("set" + grid[-1])
                if in_scope(t, g, mh)]
    if grid == "kurucz":
        from ldhybrid.grids import exoticld
        out = []
        for p in sorted((exoticld.ROOT / "kurucz").glob("MH*/teff*/logg*/kurucz_spectra.dat")):
            mh = float(p.parts[-4][2:])
            t = int(p.parts[-3][4:])
            g = float(p.parts[-2][4:])
            if in_scope(t, g, mh):
                out.append((t, g, mh))
        return sorted(out)
    raise ValueError(grid)


def load(grid, key):
    from ldhybrid.grids import exoticld, mpsatlas, phoenix
    t, g, mh = key
    if grid == "phoenix":
        return phoenix.load(t, g, mh)
    if grid in ("mps1", "mps2"):
        return mpsatlas.load("set" + grid[-1], mh, t, g)
    return exoticld.load("kurucz", mh, t, g)


def run_chunk(grid, idx, keys):
    from ldhybrid.products import model_products
    path = OUT / grid / f"chunk_{idx:04d}.npz"
    res = [model_products(load(grid, k), BANDS[grid]) for k in keys]
    arrays = {k: np.stack([r[k] for r in res]) for k in res[0]}
    arrays["keys"] = np.array([[str(x) for x in k] for k in keys])
    arrays["bands"] = np.array(BANDS[grid])
    tmp = path.with_name(path.stem + ".tmp.npz")
    np.savez_compressed(tmp, **arrays)
    tmp.replace(path)
    return idx


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("grid")
    ap.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) - 2))
    ap.add_argument("--limit", type=int, default=None)
    a = ap.parse_args()
    keys = models(a.grid)
    chunks = [keys[i:i + CHUNK] for i in range(0, len(keys), CHUNK)][: a.limit]
    (OUT / a.grid).mkdir(parents=True, exist_ok=True)
    todo = [(i, c) for i, c in enumerate(chunks) if not (OUT / a.grid / f"chunk_{i:04d}.npz").exists()]
    print(f"{a.grid}: {len(keys)} models, {len(chunks)} chunks, {len(todo)} to do, "
          f"{len(BANDS[a.grid])} bands", flush=True)
    t0 = time.time()
    with ProcessPoolExecutor(a.workers) as ex:
        futs = [ex.submit(run_chunk, a.grid, i, c) for i, c in todo]
        for n, f in enumerate(as_completed(futs), 1):
            f.result()
            el = time.time() - t0
            print(f"  {n}/{len(todo)} chunks  {el / 60:.1f} min  eta {el / n * (len(todo) - n) / 60:.0f} min",
                  flush=True)
    print("done", flush=True)


if __name__ == "__main__":
    sys.exit(main())
