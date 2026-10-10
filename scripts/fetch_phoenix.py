"""Download and reduce the in-scope PHOENIX SpecInt v3.0 grid (resumable).

    .venv/bin/python scripts/fetch_phoenix.py [--workers 4] [--limit N] [--subset-only]
"""

import argparse
import json
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

from ldhybrid.grids import phoenix
from ldhybrid.paths import REDUCED


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--subset-only", action="store_true")
    a = ap.parse_args()

    cat = phoenix.remote_catalogue()
    (REDUCED / "phoenix").mkdir(parents=True, exist_ok=True)
    with open(REDUCED / "phoenix" / "catalogue.json", "w") as f:
        json.dump(cat, f, indent=0)
    missing = set(phoenix.P2_SUBSET) - {(e["teff"], e["logg"], e["mhs"]) for e in cat}
    for m in sorted(missing):
        print("P2 subset model not on server (grid hole):", m)
    if a.subset_only:
        cat = [e for e in cat if phoenix.keep_raw(e)]
    # subset first, so the P2 study can start before the bulk finishes
    cat.sort(key=lambda e: (not phoenix.keep_raw(e), e["mhs"] != "-0.0", e["teff"], e["logg"]))
    if a.limit:
        cat = cat[: a.limit]
    print(f"{len(cat)} models ({sum(e['size'] for e in cat) / 1e9:.1f} GB raw)", flush=True)

    t0, n, fails = time.time(), 0, []
    with ThreadPoolExecutor(a.workers) as ex:
        futs = {ex.submit(phoenix.fetch_one, e): e for e in cat}
        for fut in as_completed(futs):
            e = futs[fut]
            n += 1
            try:
                status = fut.result()
            except Exception as err:
                status = f"FAIL {err}"
                fails.append(e["name"])
            if n % 25 == 0 or status.startswith("FAIL"):
                rate = n / (time.time() - t0)
                print(f"[{n}/{len(cat)}] {e['name']} {status}  "
                      f"{rate * 60:.1f} models/min", flush=True)
    print(f"finished: {n - len(fails)} ok, {len(fails)} failed", flush=True)
    for f in fails:
        print("  failed:", f)
    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    main()
