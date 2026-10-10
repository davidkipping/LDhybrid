"""Download the archive-distributed grids into data/raw/ (resumable).

    .venv/bin/python scripts/fetch_archives.py [name ...]

Names: muram, mpsatlas, stagger, kurucz_exoticld, sensitivity, plato.
With no names, everything is fetched (small files first).
"""

import json
import sys
import urllib.request

from ldhybrid.download import fetch
from ldhybrid.paths import RAW

EDMOND = "https://edmond.mpg.de/api"
ZENODO = "https://zenodo.org/api/records"


def edmond_files(doi):
    url = f"{EDMOND}/datasets/:persistentId/?persistentId=doi:{doi}"
    with urllib.request.urlopen(url, timeout=120) as r:
        d = json.load(r)
    return [(f["dataFile"]["filename"], f["dataFile"]["id"], f["dataFile"].get("filesize"))
            for f in d["data"]["latestVersion"]["files"]]


def zenodo_files(record):
    with urllib.request.urlopen(f"{ZENODO}/{record}", timeout=120) as r:
        d = json.load(r)
    return [(f["key"], f["links"]["self"], f["size"]) for f in d["files"]]


def get_edmond(doi, subdir, source, only=None):
    for name, fid, size in edmond_files(doi):
        if only and name not in only:
            continue
        print(f"  {subdir}/{name} ({(size or 0) / 1e6:.1f} MB)", flush=True)
        fetch(f"{EDMOND}/access/datafile/{fid}", RAW / subdir / name, source=source,
              expected_size=size, timeout=3600)


def get_zenodo(record, subdir, source, only=None):
    for name, url, size in zenodo_files(record):
        if only and name not in only:
            continue
        print(f"  {subdir}/{name} ({size / 1e6:.1f} MB)", flush=True)
        fetch(url, RAW / subdir / name, source=source, expected_size=size, timeout=3600)


JOBS = {
    # MURaM 3D-MHD magnetic LD library (Kostogryz et al. 2026)
    "muram": lambda: get_edmond("10.17617/3.FBTIYY", "muram", "muram-kostogryz2026"),
    # ExoTiC-LD v3.1.2 bundle: Stagger 3D (Magic et al. 2015), throughputs, Kurucz subset
    "sensitivity": lambda: get_zenodo(7874921, "exoticld", "exoticld-3.1.2",
                                      only={"Sensitivity_files.zip"}),
    "stagger": lambda: get_zenodo(7874921, "exoticld", "exoticld-3.1.2", only={"stagger.zip"}),
    "kurucz_exoticld": lambda: get_zenodo(7874921, "exoticld", "exoticld-3.1.2",
                                          only={"kurucz.zip"}),
    # PLATO consortium MARCS profiles and coefficients (Morello et al. 2022)
    "plato": lambda: get_zenodo(7339706, "plato", "plato-morello2022"),
    # MPS-ATLAS library, Set 1 and Set 2 (Kostogryz et al. 2022; Witzke et al. 2021)
    "mpsatlas": lambda: get_edmond("10.17617/3.NJ56TR", "mpsatlas", "mpsatlas-kostogryz2022"),
}


def main(names):
    for n in names or list(JOBS):
        print(f"[{n}]", flush=True)
        JOBS[n]()
    print("done", flush=True)


if __name__ == "__main__":
    main(sys.argv[1:])
