"""Fetch passband curves: SVO photometric filters, and unpack the ExoTiC-LD
HST/JWST throughputs (fetch_archives.py sensitivity must have run).

    .venv/bin/python scripts/fetch_bands.py
"""

import zipfile

from ldhybrid import bands
from ldhybrid.download import fetch, record
from ldhybrid.paths import RAW

if __name__ == "__main__":
    for name, fid in bands.SVO.items():
        print("  SVO", name, fid, flush=True)
        fetch(bands.SVO_URL.format(fid), bands.svo_path(name), source="svo-fps")
    for name, f in bands.IRAC.items():
        print("  IRAC", name, f, flush=True)
        fetch(bands.IRAC_URL.format(f), RAW / "bands" / "irac" / f, source="spitzer-ssc-irac")
    z = RAW / "exoticld" / "Sensitivity_files.zip"
    out = RAW / "bands" / "exoticld"
    out.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(z) as zf:
        for m in zf.namelist():
            if m.endswith(".csv") and not m.startswith("__MACOSX"):
                p = out / m.split("/")[-1]
                p.write_bytes(zf.read(m))
                record("exoticld-3.1.2/Sensitivity_files.zip", str(z), p)
    for name in bands.BANDS:
        lam, S = bands.response(name)
        print(f"  {name:28s} {lam[S > 0.01 * S.max()].min() / 1e4:6.3f}–"
              f"{lam[S > 0.01 * S.max()].max() / 1e4:6.3f} µm  ({len(lam)} pts)")
