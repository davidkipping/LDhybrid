"""Data locations. Everything under data/ is git-ignored; override the root
with the LDHYBRID_DATA environment variable."""

import os
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
DATA = Path(os.environ.get("LDHYBRID_DATA", REPO / "data"))
RAW = DATA / "raw"            # files exactly as downloaded
REDUCED = DATA / "reduced"    # our binned / converted products
PROVENANCE = DATA / "provenance.jsonl"
