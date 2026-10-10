"""P5: CDS/VizieR release — fixed-width ASCII tables and a byte-by-byte ReadMe,
in the style of Claret & Bloemen (2011, J/A+A/529/A75): one table per grid ×
law, one row per (model, passband, method).

    .venv/bin/python scripts/p5_cds.py

Output: data/release/cds/{ReadMe, <grid>_<law>.dat.gz}
"""

import gzip
import textwrap
from pathlib import Path

import h5py
import numpy as np

from ldhybrid.paths import DATA
from ldhybrid.products import LAWS, METHODS

SRC = DATA / "products"
OUT = DATA / "release" / "cds"
GRIDS = {"phoenix": "PHOENIX-ACES-AGSS-COND SpecInt v3.0 (Husser et al. 2013), spherical",
         "mps1": "MPS-ATLAS Set 1 (Kostogryz et al. 2022), plane-parallel",
         "mps2": "MPS-ATLAS Set 2 (Kostogryz et al. 2022), plane-parallel",
         "kurucz": "ATLAS9 (Kurucz; ExoTiC-LD subset), plane-parallel, JWST MIRI LRS only"}
COEF = {"hybrid2": ["w1", "w2"], "hybrid4": ["w1", "w2", "w3", "w4"],
        "hybrid5": ["w1", "w2", "w3", "w4", "w5"], "quadratic": ["u1", "u2"],
        "claret4": ["a1", "a2", "a3", "a4"], "power2": ["c", "alpha"]}
LAWDEF = {
    "hybrid2": "I/I(1) = 1 - w1(1-mu^2) - w2 P(mu; 0.36)",
    "hybrid4": "I/I(1) = 1 - w1(1-mu^4) - w2(1-mu^2)^2 - w3 P(mu; e0) - w4 P(mu; e1), "
               "e = 0.011762, 0.30408",
    "hybrid5": "I/I(1) = 1 - w1(1-mu^4) - w2(1-mu^2)^2 - sum_k w(3+k) P(mu; ek), "
               "e = 0.0015966, 0.041275, 0.38743",
    "quadratic": "I/I(1) = 1 - u1(1-mu) - u2(1-mu)^2",
    "claret4": "I/I(1) = 1 - sum_k ak (1-mu^(k/2)), k = 1..4",
    "power2": "I/I(1) = 1 - c(1-mu^alpha)",
}
BANDW = 25
# (label, format, width, units, explanation)
FIXED = [("Teff", "I5", 5, "K", "Effective temperature"),
         ("logg", "F4.2", 4, "[cm/s2]", "Surface gravity"),
         ("[M/H]", "F5.2", 5, "[-]", "Metallicity"),
         ("Band", f"A{BANDW}", BANDW, "---", "Passband (1)"),
         ("Meth", "A2", 2, "---", "[r/mu] Fitting method (2)"),
         ("mucri", "F8.6", 8, "---", "Limb cut mu_cri = sqrt(1 - r_c^2) (3)")]
TAIL = [("sigma", "E9.3", 9, "---", "Weighted rms of the intensity residual (4)"),
        ("dImax", "E9.3", 9, "---", "Max |Delta I/I(1)| over the tabulated mu range"),
        ("dFlux", "E10.3", 10, "---", "Relative error in the integrated flux (5)"),
        ("LCppm", "F8.2", 8, "ppm", "Worst light-curve error (6)"),
        ("Phys", "I1", 1, "---", "[0/1] 1 if positive and monotonic"),
        ("Bound", "I1", 1, "---", "[0/1] 1 if a physical constraint binds")]


def columns(law):
    coef = [(c, "F10.6", 10, "---", f"Coefficient {c}") for c in COEF[law]]
    return FIXED + coef + TAIL


def fmt(a, f, w):
    t, spec = f[0], f[1:]
    if t == "I":
        s = np.char.mod(f"%{w}d", np.asarray(a, int))
    elif t == "F":
        s = np.char.mod(f"%{w}.{spec.split('.')[1]}f", np.asarray(a, float))
    elif t == "E":
        s = np.char.mod(f"%{w}.{spec.split('.')[1]}E", np.asarray(a, float))
    else:
        s = np.char.ljust(np.asarray(a, str), w)
    if np.any(np.char.str_len(s) != w):
        raise ValueError(f"field overflow for {f}")
    return s


def write_table(grid, law):
    h = h5py.File(SRC / f"{grid}.h5")
    bands = np.array([b.decode() for b in h["bands"][()]])
    cov = h["covered"][()]
    i, j = np.nonzero(cov)
    cols_all = []
    for m in METHODS:
        k = f"{law}|{m}"
        W = h[k + "|w"][()][i, j]
        rc = h["r_c"][()][i, j]
        vals = [h["teff"][()][i], h["logg"][()][i], h["mh"][()][i], bands[j],
                np.full(len(i), m), np.sqrt(np.clip(1 - rc ** 2, 0, None))]
        vals += [W[:, n] for n in range(W.shape[1])]
        lc = h[k + "|lc_ppm"][()][i, j]
        vals += [h[k + "|rms_I"][()][i, j], h[k + "|max_abs_I"][()][i, j],
                 h[k + "|flux_err"][()][i, j], np.where(np.isfinite(lc), np.minimum(lc, 99999.99), 99999.99),
                 h[k + "|physical"][()][i, j].astype(int), h[k + "|bound"][()][i, j].astype(int)]
        cols_all.append(vals)
    spec = columns(law)
    path = OUT / f"{grid}_{law}.dat.gz"
    n = 0
    with gzip.open(path, "wt") as f:
        for vals in cols_all:
            parts = [fmt(v, s[1], s[2]) for v, s in zip(vals, spec)]
            line = parts[0]
            for p in parts[1:]:
                line = np.char.add(np.char.add(line, " "), p)
            f.write("\n".join(line.tolist()) + "\n")
            n += len(line)
    width = sum(s[2] for s in spec) + len(spec) - 1
    return path.name.removesuffix(".gz"), width, n


PAPER = Path(__file__).resolve().parent.parent / "paper" / "ldhybrid.tex"


def _groups(text, start):
    """Consecutive {...} groups after position `start` (brace-balanced)."""
    out, i = [], start
    while True:
        while i < len(text) and text[i] in " \t\n":
            i += 1
        if i >= len(text) or text[i] != "{":
            return out
        depth, j = 0, i
        while True:
            depth += {"{": 1, "}": -1}.get(text[j], 0)
            if depth == 0:
                break
            j += 1
        out.append(text[i + 1:j])
        i = j + 1


def latex_to_text(s):
    """Plain ASCII for the CDS ReadMe (the macros used in the abstract/title)."""
    import re
    s = s.split("\\thanks")[0]
    for a, b in (("\\ppm{}", " ppm"), ("\\ppm", " ppm"), ("\\varepsilon", "epsilon"),
                 ("\\mu", "mu"), ("\\teff", "Teff"), ("\\,", " "), ("--", "-"), ("~", " "),
                 ("^{-2}", "^-2")):
        s = s.replace(a, b)
    s = s.replace("$", "")
    s = re.sub(r"\\[a-zA-Z]+\s*", "", s).replace("{", "").replace("}", "")
    return re.sub(r"\s+", " ", s).strip()


def paper_title_abstract():
    tex = PAPER.read_text()
    title = latex_to_text(_groups(tex, tex.index("\\title") + len("\\title"))[0])
    parts = _groups(tex, tex.index("\\abstract") + len("\\abstract"))
    labels = ("Context", "Aims", "Methods", "Results", "Conclusions")
    abstract = [f"{lab}. {latex_to_text(p)}" for lab, p in zip(labels, parts)]
    return title, abstract


def readme(entries):
    L = []
    a = L.append
    title, abstract = paper_title_abstract()
    a("J/A+A/???/???   Limb darkening: rational hybrid laws      (<authors>, 2026)")
    a("=" * 80)
    for ln in textwrap.wrap(title, 80):
        a(ln)
    a("    <authors>")
    a("    <Astron. Astrophys. ???, ??? (2026)>")
    a("=" * 80)
    a("ADC_Keywords: Models, atmosphere ; Stars, late-type")
    a("Keywords: stars: atmospheres - planets and satellites: general -")
    a("          techniques: photometric - methods: numerical")
    a("")
    a("Abstract:")
    for para in abstract:
        for ln in textwrap.wrap(para, 76, break_on_hyphens=False):
            a("    " + ln)
    a("")
    a("Description:")
    a("    Coefficients fitted to photon-weighted specific-intensity profiles of")
    a("    spherical PHOENIX and plane-parallel MPS-ATLAS (Sets 1 and 2) and ATLAS9")
    a("    models, for passbands used in transit observations (Kepler, TESS,")
    a("    CHEOPS, NGTS, Spitzer/IRAC, Johnson-Cousins, Sloan, 2MASS, HST WFC3 and")
    a("    STIS, JWST NIRISS, NIRSpec, NIRCam and MIRI). The hybrid laws give")
    a("    elementary transit light curves for planets of any elliptical outline.")
    a("    Fits are least squares uniform in projected radius r (Meth = r) and")
    a("    uniform in mu (Meth = mu); the hybrid laws are constrained to their")
    a("    physical regions, the quadratic and 4-term laws to positive monotonic")
    a("    profiles; the power-2 law is unconstrained.")
    a("")
    a("    The laws (P(mu; e) = [(mu^2+e)^-2 - (1+e)^-2] / [e^-2 - (1+e)^-2]):")
    for law in LAWS:
        for k, ln in enumerate(textwrap.wrap(LAWDEF[law], 64)):
            a(f"      {law if k == 0 else '':9s} {ln}")
    a("")
    a("File Summary:")
    a("-" * 80)
    a(" FileName                  Lrecl   Records   Explanations")
    a("-" * 80)
    a(" ReadMe                       80         .   This file")
    for name, width, n, grid, law in entries:
        a(f" {name:25s} {width:5d} {n:9d}   {law}, {grid}")
    a("-" * 80)
    a("")
    a("Description of the models:")
    for g, d in GRIDS.items():
        for k, ln in enumerate(textwrap.wrap(d, 66)):
            a(f"    {g if k == 0 else '':7s} {ln}")
    a("")
    for law in LAWS:
        a(f"Byte-by-byte Description of file: *_{law}.dat")
        a("-" * 80)
        a("   Bytes Format Units   Label  Explanations")
        a("-" * 80)
        pos = 1
        for lab, f, w, u, e in columns(law):
            a(f"{pos:4d}-{pos + w - 1:3d}  {f:6s} {u:7s} {lab:6s} {e}")
            pos += w + 1
        a("-" * 80)
        a("")
    a("Note (1): passband transmission curves: SVO Filter Profile Service")
    a("    (photometric), Spitzer Science Center (IRAC, subarray unless _full),")
    a("    ExoTiC-LD v3.1.2 (HST and JWST modes). Integration is photon-weighted.")
    a("Note (2): r = least squares uniform in projected radius (recommended);")
    a("    mu = uniform in mu (as Claret's LSM).")
    a("Note (3): spherical models only: the stellar edge r_c is at the node of")
    a("    maximum |dI/dr| (method M1 of Claret & Southworth 2023), mu rescaled;")
    a("    0 for plane-parallel models.")
    a("Note (4): in the method's own metric (uniform in r or mu).")
    a("Note (5): (integral of I mu dmu, law) / (same, model) - 1.")
    a("Note (6): max |Delta F| between the light curve of the fitted law and the")
    a("    exact light curve of the model profile, planet radius 0.1, impact")
    a("    parameters 0, 0.3, 0.6, 0.9.")
    a("")
    a("=" * 80)
    a("(End)")
    return "\n".join(L) + "\n"


if __name__ == "__main__":
    import sys
    OUT.mkdir(parents=True, exist_ok=True)
    entries = []
    if "--readme-only" in sys.argv:
        for grid in GRIDS:
            for law in LAWS:
                with gzip.open(OUT / f"{grid}_{law}.dat.gz", "rt") as f:
                    first = f.readline().rstrip("\n")
                    n = 1 + sum(1 for _ in f)
                entries.append((f"{grid}_{law}.dat", len(first), n, grid, law))
        (OUT / "ReadMe").write_text(readme(entries))
        print("ReadMe rewritten")
        sys.exit()
    for grid in GRIDS:
        for law in LAWS:
            name, width, n = write_table(grid, law)
            entries.append((name, width, n, grid, law))
            print(f"{name:28s} Lrecl {width}  {n} records", flush=True)
    (OUT / "ReadMe").write_text(readme(entries))
    print("ReadMe written")
