"""Recompute the coefficient spectrum of ONE law in every PHOENIX
coefficient-spectrum file, in place (atomic), leaving the other laws as they are.

    .venv/bin/python scripts/p3_spectra_update.py hybrid4 [--workers N]
"""

import argparse
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np

from ldhybrid import spectra
from ldhybrid.grids import phoenix
from ldhybrid.paths import DATA

OUT = DATA / "products" / "spectra" / "phoenix"


def update(path, law):
    z = dict(np.load(path))
    mhs = "-0.0" if float(z["mh"]) == 0 else f"{float(z['mh']):+.1f}"
    cs = spectra.model_spectrum(phoenix.load(int(z["teff"]), float(z["logg"]), mhs))
    if not np.isclose(cs.r_c, float(z["r_c"]), rtol=0, atol=1e-12):
        raise RuntimeError(f"{path.name}: r_c changed ({z['r_c']} -> {cs.r_c})")
    z[f"a_{law}"] = cs.a[law].astype(np.float32)
    tmp = path.with_name(path.stem + ".tmp.npz")
    np.savez_compressed(tmp, **z)
    tmp.replace(path)
    return path.name


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("law")
    ap.add_argument("--workers", type=int, default=6)
    a = ap.parse_args()
    mark = OUT / f".{a.law}_updated"
    done = set(mark.read_text().split()) if mark.exists() else set()
    todo = [p for p in sorted(OUT.glob("*.npz")) if ".tmp" not in p.name and p.name not in done]
    print(f"{len(todo)} spectra to update", flush=True)
    with ProcessPoolExecutor(a.workers) as ex, open(mark, "a") as fm:
        futs = [ex.submit(update, p, a.law) for p in todo]
        for n, f in enumerate(as_completed(futs), 1):
            fm.write(f.result() + "\n")
            fm.flush()
            if n % 200 == 0:
                print(f"  {n}/{len(todo)}", flush=True)
    print("ALL DONE", flush=True)


if __name__ == "__main__":
    sys.exit(main())
