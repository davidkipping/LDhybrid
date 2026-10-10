"""Mirror the reference coefficient catalogues from CDS into data/raw/vizier/.

    .venv/bin/python scripts/fetch_tables.py
"""

import re
import urllib.request

from ldhybrid.download import fetch
from ldhybrid.paths import RAW

CDS = "https://cdsarc.cds.unistra.fr/ftp/J/A+A"
CATALOGUES = {
    "529/A75": "Claret & Bloemen 2011: ATLAS/PHOENIX, Kepler/CoRoT/Spitzer/photometric",
    "600/A30": "Claret 2017: TESS, ATLAS + PHOENIX (quasi-spherical, r-method)",
    "618/A20": "Claret 2018: TESS/Kepler/CoRoT/MOST, PHOENIX full spherical method",
    "664/A128": "Claret & Southworth 2022: power-2, ATLAS",
    "674/A63": "Claret & Southworth 2023: power-2 + 4-term, PHOENIX-COND M1/M2",
    "699/A97": "Claret, Hauschildt & Torres 2025: JWST, PHOENIX high resolution",
    "666/A60": "Kostogryz et al. 2022: MPS-ATLAS CLVs + coefficients (Kepler/TESS/CHEOPS/PLATO)",
}


def listing(url):
    with urllib.request.urlopen(url, timeout=120) as r:
        html = r.read().decode()
    names = re.findall(r'href="([^"/?:]+)"', html)
    return [n for n in names if not n.startswith("private_") and n != "v00"]


def mirror(cat, sub=""):
    url = f"{CDS}/{cat}/{sub}"
    for name in listing(url):
        if name == "files":
            # 666/A60/files: atmosphere structures (GBs, very slow from CDS) and
            # python-format copies of the .dat tables; not needed
            continue
        dest = RAW / "vizier" / cat.replace("/", "_") / sub / name
        print("  ", dest.relative_to(RAW), flush=True)
        fetch(url + name, dest, source=f"vizier J/A+A/{cat}")


if __name__ == "__main__":
    for cat, what in CATALOGUES.items():
        print(f"J/A+A/{cat}: {what}", flush=True)
        mirror(cat)
