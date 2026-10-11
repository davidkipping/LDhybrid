"""P5: CDS/VizieR release — fixed-width ASCII tables and a byte-by-byte ReadMe,
in the style of Claret & Bloemen (2011, J/A+A/529/A75): one table per grid ×
law, one row per (model, passband, method).

    .venv/bin/python scripts/p5_cds.py

Output: data/release/cds/{ReadMe, <grid>_<law>.dat.gz, corrections.dat,
offsets.dat, systematics.dat}
"""

import csv
import gzip
import textwrap
from pathlib import Path

import h5py
import numpy as np

from ldhybrid import v1
from ldhybrid.paths import DATA
from ldhybrid.products import LAWS, METHODS
from squishierplanet.laws import hybrid2_vertices

SRC = DATA / "products"
OUT = DATA / "release" / "cds"
GRIDS = {"phoenix": "PHOENIX-ACES-AGSS-COND SpecInt v3.0 (Husser et al. 2013), spherical",
         "mps1": "MPS-ATLAS Set 1 (Kostogryz et al. 2022), plane-parallel",
         "mps2": "MPS-ATLAS Set 2 (Kostogryz et al. 2022), plane-parallel"}
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


POLES = {"hybrid2": (v1.HYBRID2_EPS,), "hybrid4": v1.HYBRID4_EPS, "hybrid5": v1.HYBRID5_EPS}


def columns(law):
    coef = [(c, "F10.6", 10, "---", f"Coefficient {c}") for c in COEF[law]]
    extra = []
    if law.startswith("hybrid"):
        extra += [(f"q{j + 1}", "F9.6", 9, "---", f"Sampling coordinate q{j + 1} (7)")
                  for j in range(len(COEF[law]))]
        extra += [(c, "E15.8", 15, "---", f"Expanded form, {c} (8)") for c in ("c0", "c2", "c4")]
        extra += [(f"d{k + 1}", "E15.8", 15, "---", f"Expanded form, d{k + 1} (8)")
                  for k in range(len(POLES[law]))]
    return FIXED + coef + extra + TAIL


def to_q(law, W):
    """Sampling coordinates q (vectorised ldhybrid.qspace.to_q)."""
    W = np.asarray(W, float)
    if law == "hybrid2":
        Vc, Vl = hybrid2_vertices(v1.HYBRID2_EPS)
        st = np.linalg.solve(np.column_stack([Vc, Vl - Vc]), W.T).T
        s_, t_ = st[:, 0], st[:, 1]
        return np.column_stack([s_ * s_, np.where(s_ != 0, t_ / np.where(s_ != 0, s_, 1), 0.0)])
    n = W.shape[1]
    Q = np.empty_like(W)
    stick = np.ones(len(W))
    for j in range(n):
        f = np.where(stick > 0, 1.0 - W[:, j] / np.where(stick > 0, stick, 1), 1.0)
        Q[:, j] = np.where(f >= 0, np.abs(f) ** (n - j), f)
        stick = np.where(f > 0, stick * f, 0.0)
    return Q


def expanded(law, W):
    """I/I(1) = c0 + c2 mu^2 + c4 mu^4 + sum_k d_k (mu^2 + e_k)^-2 (linear in w)."""
    W = np.asarray(W, float)
    eps = np.array(POLES[law])
    N = eps ** -2 - (1 + eps) ** -2
    if law == "hybrid2":
        wp = W[:, 1:]
        c0 = 1 - W[:, 0] + wp @ ((1 + eps) ** -2 / N)
        c2, c4 = W[:, 0], np.zeros(len(W))
    else:
        wp = W[:, 2:]
        c0 = 1 - W[:, 0] - W[:, 1] + wp @ ((1 + eps) ** -2 / N)
        c2, c4 = 2 * W[:, 1], W[:, 0] - W[:, 1]
    d = -wp / N
    return np.column_stack([c0, c2, c4, d])


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
        if law.startswith("hybrid"):
            Q, X = to_q(law, W), expanded(law, W)
            vals += [Q[:, n] for n in range(Q.shape[1])] + [X[:, n] for n in range(X.shape[1])]
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


# ------------------------------------------------------------ auxiliary tables
def _fw(v, f, w):
    """One fixed-width field; None or '' -> blanks (CDS null)."""
    if v is None or (isinstance(v, str) and v.strip() == "") or (isinstance(v, float) and np.isnan(v)):
        return " " * w
    out = fmt(np.array([v]), f, w)[0]
    return out


def write_aux(name, spec, rows):
    lines = [" ".join(_fw(r.get(lab), f, w) for lab, f, w, u, e in spec) for r in rows]
    (OUT / name).write_text("\n".join(lines) + "\n")
    return name, sum(s[2] for s in spec) + len(spec) - 1, len(lines)


CORR_SPEC = [("Source", "A7", 7, "---", "Simulation code (MURaM or Stagger)"),
             ("Kind", "A8", 8, "---", "[magnetic/3d] Correction kind"),
             ("State", "A5", 5, "---", "[ssd/B100G/B200G/B300G/hydro] Simulation state (9)"),
             ("Star", "A9", 9, "---", "? MURaM simulation name"),
             ("Teff", "I5", 5, "K", "Effective temperature of the simulation"),
             ("logg", "F5.3", 5, "[cm/s2]", "Surface gravity of the simulation"),
             ("[M/H]", "F5.2", 5, "[-]", "Metallicity of the simulation"),
             ("Band", f"A{BANDW}", BANDW, "---", "Passband"),
             ("Law", "A7", 7, "---", "Hybrid law"),
             ("dw1", "F9.6", 9, "---", "Weight shift Delta w1 (10)"),
             ("dw2", "F9.6", 9, "---", "Weight shift Delta w2"),
             ("dw3", "F9.6", 9, "---", "? Weight shift Delta w3"),
             ("dw4", "F9.6", 9, "---", "? Weight shift Delta w4"),
             ("dw5", "F9.6", 9, "---", "? Weight shift Delta w5"),
             ("dLC", "F7.2", 7, "ppm", "Light-curve change the shift implies (6)"),
             ("Ext", "I1", 1, "---", "[0/1] 1 if the ratio is held constant below mu=0.1"),
             ("Base", "A22", 22, "---", "1D base profile I_base (10)")]


def corrections_rows():
    out = []
    for r in csv.DictReader(open(SRC / "corrections.csv")):
        fl = r["flags"]
        if r["source"] == "Stagger":
            base = "ATLAS9 same params"
        elif r["kind"] == "magnetic":
            base = "MPS2 nearest point"
        else:
            b = fl.split("1D=MPS2(")[-1].rstrip(")").replace(",", "")
            base = "MPS2 " + b if "1D=MPS2" in fl else "MPS2"
        out.append({"Source": r["source"], "Kind": r["kind"], "State": r["state"],
                    "Star": r["star"] or None, "Teff": int(float(r["teff"])),
                    "logg": float(r["logg"]), "[M/H]": float(r["mh"]), "Band": r["band"],
                    "Law": r["law"], **{f"dw{k}": float(r[f"dw{k}"]) if r.get(f"dw{k}") else None
                                        for k in range(1, 6)},
                    "dLC": float(r["dlc_ppm"]), "Ext": int("extrapolated_below_0.1" in fl),
                    "Base": base})
    return out


OFF_SPEC = [("Regime", "A15", 15, "---", "Passband or wavelength range of the measurements"),
            ("Bands", "A26", 26, "---", "Passbands or modes the offset applies to (11)"),
            ("Grid", "A7", 7, "---", "Model grid"),
            ("Corr", "A11", 11, "---", "Correction state the offset applies to (11)"),
            ("Route", "A7", 7, "---", "[hybrid5/profile] What the model value is (12)"),
            ("Stat", "A5", 5, "---", "[h1p/h2p/ell/delta] Statistic (13)"),
            ("Offset", "F7.4", 7, "---", "Mean offset, observed - model"),
            ("e_Offset", "F6.4", 6, "---", "? Error of the mean offset"),
            ("Scatter", "F6.4", 6, "---", "Star-to-star scatter (prior width for one star)"),
            ("rho", "F6.3", 6, "---", "? Correlation of the h1p and h2p star offsets"),
            ("N", "I2", 2, "---", "Number of stars"),
            ("Ref", "A30", 30, "---", "Source of the measurements")]


def offsets_rows():
    out = []
    for r in csv.DictReader(l for l in open(SRC / "empirical_offsets.csv") if not l.startswith("#")):
        f = lambda k: float(r[k]) if r.get(k, "").strip() else None
        out.append({"Regime": r["regime"], "Bands": r["bands"], "Grid": r["grid"],
                    "Corr": r["correction"],
                    "Route": r["route"], "Stat": r["statistic"], "Offset": f("offset"),
                    "e_Offset": f("offset_err"), "Scatter": f("scatter"), "rho": f("rho_h1h2"),
                    "N": int(r["N"]), "Ref": r["source"]})
    return out


SYS_SPEC = [("Teff", "I5", 5, "K", "Effective temperature"),
            ("logg", "F4.2", 4, "[cm/s2]", "Surface gravity"),
            ("[M/H]", "F5.2", 5, "[-]", "Metallicity"),
            ("Band", f"A{BANDW}", BANDW, "---", "Passband"),
            ("Law", "A7", 7, "---", "Hybrid law")] + \
           [(f"s_w{k}", "F8.5", 8, "---", ("? " if k > 2 else "") +
             f"s.d. of w{k} over PHOENIX, MPS-ATLAS Set 1 and Set 2") for k in range(1, 6)] + \
           [("LCsys", "F8.2", 8, "ppm", "Largest pairwise light-curve difference (14)")]


def systematics_rows():
    out = []
    with gzip.open(SRC / "systematics.csv.gz", "rt") as fh:
        for r in csv.DictReader(fh):
            out.append({"Teff": int(float(r["teff"])), "logg": float(r["logg"]),
                        "[M/H]": float(r["mh"]), "Band": r["band"], "Law": r["law"],
                        **{f"s_w{k}": float(r[f"sd_w{k}"]) if r[f"sd_w{k}"] else None
                           for k in range(1, 6)},
                        "LCsys": float(r["max_pairwise_lc_ppm"])})
    return out


AUX = [("corrections.dat", CORR_SPEC, corrections_rows,
        "Delta w corrections (3D and magnetic) for the hybrid laws"),
       ("offsets.dat", OFF_SPEC, offsets_rows, "Empirical offsets (Kepler, TESS, JWST)"),
       ("systematics.dat", SYS_SPEC, systematics_rows,
        "Grid-to-grid systematic at the shared grid points")]


def constraint_stats():
    """Fractions of rows with an active constraint (bound) and, for power-2,
    unphysical, per grid, law and method."""
    out = {}
    for g in GRIDS:
        with h5py.File(SRC / f"{g}.h5") as h:
            cov = h["covered"][()]
            for law in LAWS:
                for m in METHODS:
                    k = f"{law}|{m}"
                    out[(g, law, m)] = (float(h[k + "|bound"][()][cov].mean()),
                                        float(1 - h[k + "|physical"][()][cov].mean()))
    return out


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
                 ("\\sigma", " sigma"), ("\\ell", "ell"), ("\\Delta", "Delta "),
                 ("\\times", " x "), ("\\lesssim", "<~"), ("\\gtrsim", ">~"),
                 ("\\mu", "mu"), ("\\star", "*"), ("_{\\rm p}", "_p"), ("\\teff", "Teff"), ("\\,", " "), ("--", "-"), ("~", " "),
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


def readme(entries, aux_entries):
    st = constraint_stats()
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
    a("    spherical PHOENIX and plane-parallel MPS-ATLAS (Sets 1 and 2) models,")
    a("    for passbands used in transit observations (Kepler, TESS, CHEOPS, NGTS,")
    a("    Spitzer/IRAC, Johnson-Cousins, Sloan, 2MASS, HST WFC3 and STIS, JWST")
    a("    NIRISS, NIRSpec, NIRCam and MIRI). The hybrid laws give elementary")
    a("    transit light curves for planets of any elliptical outline.")
    a("    Fits are least squares uniform in projected radius r (Meth = r) and")
    a("    uniform in mu (Meth = mu); the hybrid laws are constrained to their")
    a("    physical regions, the quadratic and 4-term laws to positive monotonic")
    a("    profiles; the power-2 law is unconstrained.")
    a("    The r-weighted coefficients (Meth = r) are optimised for transit light")
    a("    curves; they fit the extreme limb less well. For uses that depend on")
    a("    the limb intensities themselves (eclipsing binaries, interferometry,")
    a("    microlensing) use the uniform-mu fits (Meth = mu).")
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
    for name, width, n, expl in aux_entries:
        a(f" {name:25s} {width:5d} {n:9d}   {expl}")
    a("-" * 80)
    a("")
    a("Description of the models:")
    for g, d in GRIDS.items():
        for k, ln in enumerate(textwrap.wrap(d, 66)):
            a(f"    {g if k == 0 else '':7s} {ln}")
    a("")
    for g in GRIDS:
        with h5py.File(SRC / f"{g}.h5") as h:
            bl = [b.decode() for b in h["bands"][()]]
        a(f"Passbands of {g} ({len(bl)}):")
        for ln in textwrap.wrap(", ".join(bl), 74):
            a("    " + ln)
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
    for name, spec, _, _ in AUX:
        a(f"Byte-by-byte Description of file: {name}")
        a("-" * 80)
        a("   Bytes Format Units   Label    Explanations")
        a("-" * 80)
        pos = 1
        for lab, f, w, u, e in spec:
            a(f"{pos:4d}-{pos + w - 1:3d}  {f:6s} {u:7s} {lab:8s} {e}")
            pos += w + 1
        a("-" * 80)
        a("")
    a("Note (1): passband transmission curves: SVO Filter Profile Service")
    a("    (photometric), Spitzer Science Center (IRAC, subarray unless _full),")
    a("    ExoTiC-LD v3.1.2 (HST and JWST modes). Integration is photon-weighted.")
    a("    JWST_NIRSpec_G395H_NRS1/NRS2: the G395H curve restricted to the NRS1")
    a("    (2.87-3.72 um) and NRS2 (3.82-5.18 um) detector ranges (JDox, BOTS).")
    a("    PHOENIX spectra stop at 5.5 um, so MIRI LRS is given for MPS-ATLAS only.")
    a("Note (2): r = least squares uniform in projected radius (recommended for")
    a("    transits); mu = uniform in mu (as Claret's LSM; for limb-sensitive uses).")
    a("    Rows with an active constraint (Bound = 1), and unphysical power-2 rows")
    a("    (Phys = 0; c > 1, i.e. negative limb intensity), per grid and method:")
    for g in GRIDS:
        for m in METHODS:
            parts = [f"{law} {100 * st[(g, law, m)][0]:.1f}%" for law in LAWS if law != "power2"]
            a(f"      {g:7s} {m:2s}: bound " + ", ".join(parts[:3]) + ",")
            a(f"      {'':7s} {'':2s}  " + ", ".join(parts[3:]) +
              f"; power-2 unphysical {100 * st[(g, 'power2', m)][1]:.1f}%")
    a("    A binding constraint breaks the linearity of the fit in the data; power-2")
    a("    rows with c > 1 are kept as fitted (as in published tables) and flagged.")
    a("Note (3): spherical models only: the stellar edge r_c is at the node of")
    a("    maximum |dI/dr| (method M1 of Claret & Southworth 2023), mu rescaled;")
    a("    0 for plane-parallel models.")
    a("Note (4): in the method's own metric (uniform in r or mu).")
    a("Note (5): (integral of I mu dmu, law) / (same, model) - 1.")
    a("Note (6): max |Delta F| between the light curve of the fitted law and the")
    a("    exact light curve of the model profile, planet radius 0.1, impact")
    a("    parameters 0, 0.3, 0.6, 0.9.")
    a("Note (7): hybrid laws only. Coordinates on the unit hypercube that map to the")
    a("    physical region with a uniform prior: hybrid2 w = sqrt(q1)[(1-q2)Vc + q2 Vl]")
    a("    on the exact triangle with vertices O, Vc = (-0.124615, 1.124615) and")
    a("    Vl = (1.201039, -0.201039); hybrid4/5 stick-breaking of a flat Dirichlet,")
    a("    w_j = (1 - q_j^(1/(n+1-j))) prod_(i<j) q_i^(1/(n+1-i)).")
    a("Note (8): hybrid laws only. The same law in expanded form,")
    a("    I/I(1) = c0 + c2 mu^2 + c4 mu^4 + sum_k d_k (mu^2 + e_k)^-2, with the")
    a("    pole positions e_k of the law (c4 = 0 for hybrid2).")
    a("Note (9): MURaM (Voegler et al. 2005; Kostogryz et al. 2026): ssd =")
    a("    small-scale dynamo, B100G/B200G/B300G = vertical fields, hydro = field-")
    a("    free 3D; Stagger (Magic et al. 2015) 3D, as distributed with ExoTiC-LD.")
    a("Note (10): Delta w = w(I_base C) - w(I_base), added to the tabulated weights,")
    a("    with C the centre-to-limb ratio of the corrected to the base profile and")
    a("    I_base the MPS-ATLAS Set 2 band profile at the nearest grid point.")
    a("    Magnetic: C = I_B / I_hydro (MURaM with and without field). MURaM 3D:")
    a("    C = I_hydro / I_MPS2 (MPS-ATLAS also synthesised the MURaM spectra; grid")
    a("    point within 150 K, given in Base). Stagger 3D: Delta w = w(Stagger) -")
    a("    w(ATLAS9) at identical parameters (ATLAS9: Castelli & Kurucz 2003, as in")
    a("    ExoTiC-LD); adding it to MPS-ATLAS weights assumes the shift transfers")
    a("    between 1D grids. MURaM has no intensities below mu = 0.1: the ratio is")
    a("    held constant there (Ext = 1), so the shifts of the inner-pole weights")
    a("    are set by that assumption. MURaM stars: Teff = 3196, 3518, 3707, 4543,")
    a("    5240, 5778, 6807 K at [M/H] = 0 (log g 5.0-4.0), and 5750/5780 K at")
    a("    [M/H] = -1.0/+0.5. The interpolation code interpolates the [M/H] = 0")
    a("    shifts linearly in Teff only (clamped outside 3196-6807 K) and applies")
    a("    them at any log g and [M/H]: an extrapolation away from the simulated")
    a("    gravities and metallicity.")
    a("Note (11): the offset is measured against predictions in this correction")
    a("    state (none, or MURaM_B100G = Delta w of B100G added). Apply an offset only")
    a("    to a prior built in the same state, otherwise the correction is counted")
    a("    twice. Kepler and TESS rows apply to those passbands; JWST rows were")
    a("    measured over 0.86-5.5 um (NIRISS SOSS and NIRSpec PRISM bins).")
    a("Note (12): hybrid5 = the statistic of the tabulated hybrid5 fit (use these")
    a("    to condition hybrid-law weights, since the statistic is then a linear")
    a("    function of the weights); profile = the statistic of the model profile")
    a("    itself (Sing et al.'s comparison).")
    a("Note (13): h1p = I(2/3)/I(1), h2p = h1p - I(1/3)/I(1) (Maxted 2023); ell =")
    a("    1 - (c1 + c2) and delta = c2/4 of a quadratic law fitted uniformly in mu")
    a("    for mu >= 0.2 (Sing et al. 2026). rho: correlation of the per-star h1p and")
    a("    h2p offsets, for a joint (2D) Gaussian update.")
    a("Note (14): at the 1020 (Teff, log g, [M/H]) points shared by PHOENIX and")
    a("    both MPS-ATLAS sets (Teff >= 3500 K): the largest pairwise difference of")
    a("    the three grids' light curves (planet radius 0.1, as Note 6).")
    a("")
    a("=" * 80)
    a("(End)")
    return "\n".join(L) + "\n"


if __name__ == "__main__":
    import sys
    OUT.mkdir(parents=True, exist_ok=True)
    for stale in OUT.glob("kurucz_*.dat.gz"):          # ATLAS9 dropped from the release
        stale.unlink()
    entries = []
    aux_entries = []
    if "--readme-only" in sys.argv:
        for grid in GRIDS:
            for law in LAWS:
                with gzip.open(OUT / f"{grid}_{law}.dat.gz", "rt") as f:
                    first = f.readline().rstrip("\n")
                    n = 1 + sum(1 for _ in f)
                entries.append((f"{grid}_{law}.dat", len(first), n, grid, law))
        for name, spec, _, expl in AUX:
            lines = (OUT / name).read_text().splitlines()
            aux_entries.append((name, sum(x[2] for x in spec) + len(spec) - 1, len(lines), expl))
        (OUT / "ReadMe").write_text(readme(entries, aux_entries))
        print("ReadMe rewritten")
        sys.exit()
    for name, spec, rows, expl in AUX:
        nm, width, n = write_aux(name, spec, rows())
        aux_entries.append((nm, width, n, expl))
        print(f"{nm:28s} Lrecl {width}  {n} records", flush=True)
    if "--aux-only" not in sys.argv:
        for grid in GRIDS:
            for law in LAWS:
                name, width, n = write_table(grid, law)
                entries.append((name, width, n, grid, law))
                print(f"{name:28s} Lrecl {width}  {n} records", flush=True)
    else:
        for grid in GRIDS:
            for law in LAWS:
                with gzip.open(OUT / f"{grid}_{law}.dat.gz", "rt") as f:
                    first = f.readline().rstrip("\n")
                    n = 1 + sum(1 for _ in f)
                entries.append((f"{grid}_{law}.dat", len(first), n, grid, law))
    (OUT / "ReadMe").write_text(readme(entries, aux_entries))
    print("ReadMe written")
