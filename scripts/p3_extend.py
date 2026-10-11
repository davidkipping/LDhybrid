"""P3 extension: add passbands to existing production chunks, and export the
band-integrated profiles of every model × band.

    .venv/bin/python scripts/p3_extend.py phoenix mps1 mps2 [--workers N]

For each chunk in data/products/bands/<grid>/, the bands of p3_tables.BANDS
missing from it are computed with products.model_products and appended along
the band axis. In the same pass the post-limb-cut profiles (exactly as
band_profile returns them: μ′ ascending to 1, y = I/I(1)) of all bands are
written to data/products/profiles/<grid>/chunk_NNNN.npz, then gathered into

    data/products/profiles/<grid>.h5
        teff, logg, mh (n_model); bands (S);
        mu, y (n_model, n_band, n_mu_max) float32, NaN-padded; r_c (n_model, n_band)

Resumable: chunks already extended, and profile chunks already written, are
skipped.
"""

import argparse
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

import h5py
import numpy as np

from ldhybrid import bands as B
from ldhybrid.paths import DATA

sys.path.insert(0, str(DATA.parent / "scripts"))
from p3_tables import BANDS, load  # noqa: E402

OUT = DATA / "products" / "bands"
PROF = DATA / "products" / "profiles"


def parse_key(grid, k):
    t, g, mh = k
    return (int(float(t)), float(g), mh if grid == "phoenix" else float(mh))


def profiles(spec, band_names, nmu):
    from ldhybrid.profile import band_profile
    nb = len(band_names)
    mu = np.full((nb, nmu), np.nan, np.float32)
    y = np.full((nb, nmu), np.nan, np.float32)
    rc = np.full(nb, np.nan, np.float32)
    spherical = spec.meta.get("geometry") == "spherical"
    for j, b in enumerate(band_names):
        lam, S = B.response(b)
        if spec.coverage(lam, S) < 0.999:
            continue
        pr = band_profile(spec, lam, S, limb="M1") if spherical else band_profile(spec, lam, S)
        n = len(pr.mu)
        mu[j, :n], y[j, :n] = pr.mu, pr.y
        rc[j] = pr.meta.get("r_c", 1.0)
    return mu, y, rc


def run_chunk(grid, path):
    from ldhybrid.products import model_products
    z = dict(np.load(path))
    have = [str(b) for b in z["bands"]]
    missing = [b for b in BANDS[grid] if b not in have]
    ppath = PROF / grid / path.name
    if not missing and ppath.exists():
        return path.name, 0
    keys = [parse_key(grid, k) for k in z["keys"]]
    allb = have + missing
    nmu = None
    P = {"mu": [], "y": [], "r_c": []}
    new = []
    for k in keys:
        spec = load(grid, k)
        if nmu is None:
            nmu = len(spec.mu) + 1
        if missing:
            new.append(model_products(spec, missing))
        if not ppath.exists():
            m, yy, rc = profiles(spec, allb, nmu)
            P["mu"].append(m); P["y"].append(yy); P["r_c"].append(rc)
    if not ppath.exists():
        ppath.parent.mkdir(parents=True, exist_ok=True)
        tmp = ppath.with_name(ppath.stem + ".tmp.npz")
        np.savez_compressed(tmp, keys=z["keys"], bands=np.array(allb),
                            **{k: np.stack(v) for k, v in P.items()})
        tmp.replace(ppath)
    if missing:
        add = {k: np.stack([r[k] for r in new]) for k in new[0]}
        for k, v in add.items():
            z[k] = np.concatenate([z[k], v], axis=1)
        z["bands"] = np.array(allb)
        tmp = path.with_name(path.stem + ".tmp.npz")
        np.savez_compressed(tmp, **z)
        tmp.replace(path)
    return path.name, len(missing)


def gather_profiles(grid):
    files = sorted((PROF / grid).glob("chunk_*.npz"))
    parts = [dict(np.load(f)) for f in files]
    bands = parts[0]["bands"]
    assert all(list(p["bands"]) == list(bands) for p in parts)
    nmu = max(int(np.max(np.sum(np.isfinite(p["mu"]), axis=-1))) for p in parts)
    keys = np.concatenate([p["keys"] for p in parts])
    with h5py.File(PROF / f"{grid}.h5", "w") as h:
        h["teff"] = keys[:, 0].astype(float)
        h["logg"] = keys[:, 1].astype(float)
        h["mh"] = keys[:, 2].astype(float)
        h["bands"] = bands.astype("S")
        for k in ("mu", "y"):
            h.create_dataset(k, data=np.concatenate([p[k][..., :nmu] for p in parts]),
                             compression="gzip", shuffle=True)
        h["r_c"] = np.concatenate([p["r_c"] for p in parts])
        h.attrs["layout"] = ("post-limb-cut band profiles: mu ascending to 1, y = I/I(1), "
                             "NaN-padded; [model, band, node]")
    print(f"{grid}: profiles {len(keys)} models × {len(bands)} bands, n_mu ≤ {nmu}", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("grids", nargs="+")
    ap.add_argument("--workers", type=int, default=5)
    a = ap.parse_args()
    for grid in a.grids:
        files = sorted((OUT / grid).glob("chunk_*.npz"))
        files = [f for f in files if ".tmp" not in f.name]
        print(f"{grid}: {len(files)} chunks, target {len(BANDS[grid])} bands", flush=True)
        t0 = time.time()
        with ProcessPoolExecutor(a.workers) as ex:
            futs = [ex.submit(run_chunk, grid, f) for f in files]
            for n, fu in enumerate(as_completed(futs), 1):
                fu.result()
                if n % 10 == 0 or n == len(files):
                    el = time.time() - t0
                    print(f"  {n}/{len(files)}  {el / 60:.1f} min  eta "
                          f"{el / n * (len(files) - n) / 60:.0f} min", flush=True)
        gather_profiles(grid)
    print("done", flush=True)


if __name__ == "__main__":
    sys.exit(main())
