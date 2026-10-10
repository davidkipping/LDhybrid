"""Shared set-up for the P2 method studies: the representative subset, the
bands, and a parallel runner that writes JSON lines.

Subset:
  * PHOENIX v3.0 (spherical, M1 cut unless a study varies it): the 34 P2
    models with raw files kept (grids/phoenix.P2_SUBSET);
  * MPS-ATLAS Set 2 (plane-parallel): an FGK sequence at [M/H] = 0, log g 4.5,
    plus [M/H] = −1, +0.5 at three temperatures.
Bands: two photometric, two HST, four JWST — blue to 5 µm.
"""

from __future__ import annotations

import json
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np

BANDS = ["TESS", "Kepler", "HST_STIS_G430L", "HST_WFC3_G141", "JWST_NIRISS_SOSS_o1",
         "JWST_NIRSpec_G235H", "JWST_NIRSpec_G395H", "JWST_NIRSpec_PRISM"]

PHOENIX = None   # filled lazily from grids.phoenix.P2_SUBSET
MPS = ([(0.0, t, 4.5) for t in (3600, 4000, 4400, 4800, 5200, 5600, 6000, 6400, 6800, 7200)]
       + [(mh, t, 4.5) for t in (4400, 5600, 6800) for mh in (-1.0, 0.5)])


def models(which=("phoenix", "mps")):
    from ldhybrid.grids import phoenix
    out = []
    if "phoenix" in which:
        out += [("phoenix", (t, g, mh)) for t, g, mh in phoenix.P2_SUBSET]
    if "mps" in which:
        out += [("mps2", m) for m in MPS]
    return out


def load(kind, key):
    from ldhybrid.grids import mpsatlas, phoenix
    if kind == "phoenix":
        t, g, mh = key
        return phoenix.load(t, g, mh)
    if kind == "mps2":
        mh, t, g = key
        return mpsatlas.load("set2", mh, t, g)
    raise ValueError(kind)


def run_parallel(task, jobs, out_path: Path, workers: int | None = None):
    """task(job) -> list[dict]; results appended to out_path as JSON lines
    (resumable: jobs whose 'job_id' is already present are skipped)."""
    done = set()
    if out_path.exists():
        for line in open(out_path):
            done.add(json.loads(line)["job_id"])
    todo = [j for j in jobs if j["job_id"] not in done]
    print(f"{len(todo)} jobs ({len(done)} already done)", flush=True)
    workers = workers or max(1, (os.cpu_count() or 2) - 2)
    with ProcessPoolExecutor(workers) as ex, open(out_path, "a") as f:
        futs = [ex.submit(task, j) for j in todo]
        for n, fut in enumerate(as_completed(futs), 1):
            for row in fut.result():
                f.write(json.dumps(row) + "\n")
            f.flush()
            if n % 20 == 0:
                print(f"  {n}/{len(todo)}", flush=True)


def read(path: Path) -> list[dict]:
    return [json.loads(l) for l in open(path)]


def table(rows, key, laws, stat=np.median):
    return {l: float(stat([r[key] for r in rows if r["law"] == l])) for l in laws}
