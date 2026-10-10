"""Downloads with retries, checksums and a provenance log.

Every file that enters data/ is recorded as one JSON line in
data/provenance.jsonl: source, URL, local path, size, sha256, UTC time.
"""

from __future__ import annotations

import hashlib
import json
import shutil
import time
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

from .paths import PROVENANCE

USER_AGENT = "ldhybrid/0.1 (research; limb-darkening tables)"


def sha256(path: Path, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while b := f.read(chunk):
            h.update(b)
    return h.hexdigest()


def record(source: str, url: str, path: Path, **extra) -> dict:
    path = Path(path)
    entry = dict(source=source, url=url, path=str(path), size=path.stat().st_size,
                 sha256=sha256(path),
                 utc=datetime.now(timezone.utc).isoformat(timespec="seconds"),
                 **extra)
    PROVENANCE.parent.mkdir(parents=True, exist_ok=True)
    with open(PROVENANCE, "a") as f:
        f.write(json.dumps(entry) + "\n")
    return entry


def fetch(url: str, dest: Path, *, source: str, expected_size: int | None = None,
          retries: int = 4, timeout: float = 600, log: bool = True) -> Path:
    """Download url to dest (atomically, via dest.part). Skips the download if
    dest exists with the expected size. Returns dest."""
    dest = Path(dest)
    if dest.exists() and (expected_size is None or dest.stat().st_size == expected_size):
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    part = dest.with_name(dest.name + ".part")
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
            with urllib.request.urlopen(req, timeout=timeout) as r, open(part, "wb") as f:
                shutil.copyfileobj(r, f, 1 << 20)
            size = part.stat().st_size
            if expected_size is not None and size != expected_size:
                raise IOError(f"size {size} != expected {expected_size}")
            part.replace(dest)
            if log:
                record(source, url, dest)
            return dest
        except Exception as e:  # network errors, truncated transfers
            if attempt == retries - 1:
                part.unlink(missing_ok=True)
                raise RuntimeError(f"failed to fetch {url}: {e}") from e
            time.sleep(5 * 2 ** attempt)
    return dest


def list_ftp(url: str, timeout: float = 120) -> list[tuple[str, int]]:
    """(name, size) pairs from a plain FTP directory listing."""
    req = urllib.request.Request(url if url.endswith("/") else url + "/")
    with urllib.request.urlopen(req, timeout=timeout) as r:
        lines = r.read().decode().splitlines()
    out = []
    for line in lines:
        parts = line.split()
        if len(parts) >= 9 and not line.startswith("d"):
            out.append((parts[-1], int(parts[4])))
    return out
