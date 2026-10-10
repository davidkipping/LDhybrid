"""P3: coefficient spectra (hybrid2/4/5) for every PHOENIX model, R = 3000.

    .venv/bin/python scripts/p3_spectra.py

Writes data/products/spectra/phoenix/<model>.npz with lam_edges, I1, a_<law>
(float32 except the edges), r_c and the model parameters.
"""

import json

import numpy as np

from ldhybrid import spectra
from ldhybrid.grids import phoenix
from ldhybrid.paths import DATA, REDUCED

OUT = DATA / "products" / "spectra" / "phoenix"

if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    cat = json.load(open(REDUCED / "phoenix" / "catalogue.json"))
    for n, e in enumerate(sorted(cat, key=lambda e: (e["mhs"], e["teff"], e["logg"]))):
        path = OUT / e["name"].replace("SpecInt.fits", "coeffspec.npz")
        if path.exists():
            continue
        cs = spectra.model_spectrum(phoenix.load(e["teff"], e["logg"], e["mhs"]))
        np.savez_compressed(path, lam_edges=cs.lam_edges, I1=cs.I1.astype(np.float32),
                            r_c=cs.r_c, teff=e["teff"], logg=e["logg"], mh=float(e["mhs"]),
                            **{f"a_{k}": v.astype(np.float32) for k, v in cs.a.items()})
        if n % 200 == 0:
            print(n, path.name, flush=True)
    print("done")
