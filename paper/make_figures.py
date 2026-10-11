"""Figures for the LDhybrid v1.0 paper (A&A widths: 88 mm column, 180 mm page).

    .venv/bin/python paper/make_figures.py

Style: system sans (Helvetica), neutral ink and hairline chrome, focus + context.
The three hybrid laws carry the first three categorical slots of the dataviz
reference palette (validated all-pairs, light mode: blue, orange, aqua); the
control laws recede into neutral greys and are told apart by line style.
Passbands use an ordinal blue ramp ordered by wavelength (steps 250-650, each
>= 2:1 against white). Text always wears neutral ink.
"""

import json
from pathlib import Path

import h5py
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

from ldhybrid import bands, fit, laws, v1
from ldhybrid.grids import phoenix
from ldhybrid.paths import DATA
from ldhybrid.products import constraint_for
from ldhybrid.profile import band_profile

HERE = Path(__file__).parent
FIG = HERE / "figures"
FIG.mkdir(exist_ok=True)
ROOT = HERE.parent
COL, PAGE = 88 / 25.4, 180 / 25.4                      # inches

# --- tokens -----------------------------------------------------------------
INK, INK2, MUTED = "#0b0b0b", "#52514e", "#898781"
GRID, AXIS, WASH = "#e1e0d9", "#c3c2b7", "#f0efec"
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
RAMP = {250: "#86b6ef", 350: "#5598e7", 450: "#2a78d6", 550: "#1c5cab", 650: "#104281"}

STYLE = {  # law: (colour, linestyle, width, label)
    "hybrid5": (BLUE, "-", 1.6, "hybrid5"),
    "hybrid4": (ORANGE, "-", 1.6, "hybrid4"),
    "hybrid2": (AQUA, "-", 1.6, "hybrid2"),
    "claret4": (INK2, "-", 1.0, "four-term"),
    "power2": (MUTED, (0, (5, 2)), 1.0, "power-2"),
    "quadratic": (MUTED, (0, (1, 1.6)), 1.1, "quadratic"),
}
LAW_ORDER = ["hybrid5", "hybrid4", "hybrid2", "claret4", "power2", "quadratic"]

mpl.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    "mathtext.fontset": "custom", "mathtext.rm": "Helvetica", "mathtext.it": "Helvetica:italic",
    "mathtext.bf": "Helvetica:bold", "font.size": 7.5, "axes.labelsize": 7.5,
    "axes.titlesize": 7.5, "xtick.labelsize": 6.5, "ytick.labelsize": 6.5, "legend.fontsize": 6.8,
    "text.color": INK, "axes.labelcolor": INK2, "axes.edgecolor": AXIS, "axes.linewidth": 0.6,
    "xtick.color": MUTED, "ytick.color": MUTED, "xtick.labelcolor": INK2, "ytick.labelcolor": INK2,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6, "xtick.major.size": 2.5,
    "ytick.major.size": 2.5, "xtick.minor.size": 1.5, "ytick.minor.size": 1.5,
    "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True,
    "axes.grid.axis": "y", "grid.color": GRID, "grid.linewidth": 0.5,
    "axes.titlelocation": "left", "axes.titleweight": "bold", "axes.titlepad": 4,
    "legend.frameon": False, "legend.handlelength": 2.2, "legend.columnspacing": 1.2,
    "lines.solid_capstyle": "round", "savefig.dpi": 400, "savefig.bbox": "tight",
    "savefig.pad_inches": 0.02, "pdf.fonttype": 42,
})


def save(fig, name):
    for ext in ("pdf", "png"):
        fig.savefig(FIG / f"{name}.{ext}")
    plt.close(fig)
    print("wrote", name)


def law_handle(n):
    c, ls, lw, lab = STYLE[n]
    return Line2D([], [], color=c, ls=ls, lw=lw, label=lab)


def top_legend(fig, handles, y=1.0, ncol=None):
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, y),
               ncol=ncol or len(handles), handletextpad=0.5)


def clean_x_grid(a):
    a.grid(True, axis="x")
    a.grid(False, axis="y")


# ---------------------------------------------------------------------------
def fig_profiles():
    cases = [((5800, 4.5, "-0.0"), "TESS", "5800 K dwarf · TESS"),
             ((3000, 5.0, "-0.0"), "JWST_NIRSpec_G395H", "3000 K dwarf · NIRSpec G395H")]
    show = ["quadratic", "power2", "hybrid2", "claret4", "hybrid4", "hybrid5"]
    fig, ax = plt.subplots(2, 2, figsize=(PAGE, 3.5), sharex=True,
                           gridspec_kw=dict(height_ratios=[1.6, 1], hspace=0.12, wspace=0.16))
    for c, (key, band, title) in enumerate(cases):
        pr = band_profile(phoenix.load(*key), *bands.response(band), limb="M1")
        mu = np.linspace(0, 1, 800)
        y = pr.interpolant()(mu)
        top, bot = ax[0, c], ax[1, c]
        top.plot(mu, y, color=AXIS, lw=3.0, zorder=1)
        top.plot(pr.mu, pr.y, "o", ms=2.2, mfc="white", mec=INK2, mew=0.6, zorder=4)
        ins = top.inset_axes([0.47, 0.13, 0.5, 0.42])
        ins.plot(mu, y, color=AXIS, lw=3.0)
        ins.plot(pr.mu, pr.y, "o", ms=2.4, mfc="white", mec=INK2, mew=0.6, zorder=4)
        for n in show:
            if n == "power2":
                c_, a_ = fit.fit_power2(pr, "r", v1.RESAMPLE).w
                Im = laws.power2(c_, a_, mu)
            else:
                L = laws.get(n)
                Im = L.intensity(fit.fit_linear(L, pr, "r", v1.RESAMPLE, constraint_for(n)).w, mu)
            col, ls, lw, _ = STYLE[n]
            top.plot(mu, Im, color=col, ls=ls, lw=lw * 0.8, zorder=3)
            ins.plot(mu, Im, color=col, ls=ls, lw=lw, zorder=3)
            bot.plot(mu, (Im - y) * 1e3, color=col, ls=ls, lw=lw)
        ins.set_xlim(0, 0.12)
        ins.set_ylim(y[0] - 0.04, y[np.searchsorted(mu, 0.12)] + 0.03)
        ins.tick_params(labelsize=5.5)
        ins.set_title(r"limb, $\mu<0.12$", fontsize=6, fontweight="normal", color=INK2)
        for s in ins.spines.values():
            s.set_color(AXIS)
        ins.spines[["top", "right"]].set_visible(True)
        ins.grid(False)
        top.set_title(title)
        bot.axhline(0, color=INK2, lw=0.6)
        bot.set_yscale("symlog", linthresh=1, linscale=0.6)
        bot.set_ylim(-400, 400)
        bot.set_yticks([-100, -10, -1, 0, 1, 10, 100])
        bot.set_yticklabels(["−100", "−10", "−1", "0", "1", "10", "100"])
        bot.set_xlabel(r"$\mu$")
        bot.set_xlim(0, 1)
    ax[0, 0].set_ylabel(r"$I(\mu)/I(1)$")
    ax[1, 0].set_ylabel(r"law − model  ($10^{-3}$)")
    nodes = Line2D([], [], ls="", marker="o", ms=3, mfc="white", mec=INK2, mew=0.6,
                   label="PHOENIX nodes")
    interp = Line2D([], [], color=AXIS, lw=3.0, label="model (interpolated)")
    top_legend(fig, [nodes, interp] + [law_handle(n) for n in LAW_ORDER], y=0.945)
    save(fig, "fig_profiles")


def fig_quality():
    fig, ax = plt.subplots(2, 1, figsize=(COL, 4.6), sharex=True, gridspec_kw=dict(hspace=0.22))
    for a, grid, title in ((ax[0], "phoenix", "PHOENIX v3.0 (spherical)"),
                           (ax[1], "mps2", "MPS-ATLAS Set 2 (plane-parallel)")):
        h = h5py.File(DATA / "products" / f"{grid}.h5")
        T, cov = h["teff"][()], h["covered"][()]
        Tn = np.unique(T)
        for n in LAW_ORDER[::-1]:
            e = h[f"{n}|r|lc_ppm"][()].astype(float)
            e[~cov] = np.nan
            sel = [T == t for t in Tn]
            med = np.array([np.nanmedian(e[s_]) for s_ in sel])
            col, ls, lw, _ = STYLE[n]
            if n.startswith("hybrid"):
                lo = np.array([np.nanpercentile(e[s_], 16) for s_ in sel])
                hi = np.array([np.nanpercentile(e[s_], 84) for s_ in sel])
                a.fill_between(Tn, lo, hi, color=col, alpha=0.13, lw=0)
            a.plot(Tn, med, color=col, ls=ls, lw=lw)
        a.set_yscale("log")
        a.set_ylim(0.2, 300)
        a.set_xlim(2300, 7500)
        a.set_title(title)
        a.set_ylabel("light-curve error (ppm)")
        a.yaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda v, _: f"{v:g}"))
    ax[1].set_xlabel(r"$T_{\rm eff}$ (K)")
    fig.subplots_adjust(top=0.86)
    top_legend(fig, [law_handle(n) for n in LAW_ORDER], y=0.885, ncol=3)
    save(fig, "fig_quality")


def fig_basis():
    """Shapes T_j(μ) of hybrid5 (and the hybrid2 pole), one column."""
    mu = np.linspace(0, 1, 1200)
    fig, ax = plt.subplots(figsize=(COL, 2.4))
    L = laws.get("hybrid5")
    T = L.T(mu)
    labs = [r"$1-\mu^4$", r"$(1-\mu^2)^2$"] + [rf"$\Pi_{{{e:.4g}}}$" for e in v1.HYBRID5_EPS]
    cols = [INK2, MUTED, RAMP[650], RAMP[450], RAMP[250]]
    sty = ["-", (0, (4, 1.5)), "-", "-", "-"]
    for j in range(5):
        ax.plot(mu, T[:, j], color=cols[j], ls=sty[j], lw=1.4, label=labs[j])
    ax.plot(mu, laws.get("hybrid2").T(mu)[:, 1], color=AQUA, ls=(0, (1, 1.3)), lw=1.4,
            label=r"$\Pi_{0.36}$ (hybrid2)")
    ax.set_xscale("log")
    ax.set_xlim(1e-3, 1)
    ax.set_ylim(-0.02, 1.02)
    ax.set_xlabel(r"$\mu$")
    ax.set_ylabel(r"$T_j(\mu)$")
    ax.legend(loc="lower left", fontsize=6.2, handlelength=2.0)
    clean_x_grid(ax)
    save(fig, "fig_basis")


def fig_coefficients():
    bset = [("Kepler", "Kepler", 250), ("TESS", "TESS", 350),
            ("HST_WFC3_G141", "WFC3 G141", 450), ("JWST_NIRISS_SOSS_o1", "NIRISS SOSS", 550),
            ("JWST_NIRSpec_G395H", "NIRSpec G395H", 650)]
    dash = {250: (0, (1, 1.3)), 350: (0, (4, 1.5)), 450: "-", 550: (0, (6, 1.5, 1, 1.5)), 650: "-"}
    width = {250: 1.3, 350: 1.2, 450: 1.0, 550: 1.2, 650: 1.9}
    fig, ax = plt.subplots(2, 5, figsize=(PAGE, 3.3), sharex=True, sharey="col",
                           gridspec_kw=dict(hspace=0.12, wspace=0.32))
    for r, (grid, label) in enumerate((("phoenix", "PHOENIX v3.0"), ("mps2", "MPS-ATLAS Set 2"))):
        h = h5py.File(DATA / "products" / f"{grid}.h5")
        B = [b.decode() for b in h["bands"][()]]
        sel = (h["logg"][()] == 4.5) & (np.abs(h["mh"][()]) < 1e-9)
        T = h["teff"][()][sel]
        o = np.argsort(T)
        W = h["hybrid5|r|w"][()][sel][o]
        for b, lab, step in bset:
            for j in range(5):
                ax[r, j].plot(T[o], W[:, B.index(b), j], color=RAMP[step], ls=dash[step],
                              lw=width[step])
        ax[r, 0].annotate(label, xy=(-0.42, 0.5), xycoords="axes fraction", rotation=90,
                          ha="center", va="center", fontsize=7.5, fontweight="bold", color=INK)
    for j in range(5):
        ax[0, j].set_title(f"$w_{j + 1}$")
        ax[1, j].set_xlabel(r"$T_{\rm eff}$ (K)")
        ax[1, j].set_xlim(2300, 7500)
        ax[1, j].xaxis.set_major_locator(mpl.ticker.MultipleLocator(2000))
    handles = [Line2D([], [], color=RAMP[s], ls=dash[s], lw=width[s] + 0.3, label=lab)
               for _, lab, s in bset]
    top_legend(fig, handles, y=0.955)
    save(fig, "fig_coefficients")


def fig_validation():
    import csv
    mx = json.loads((ROOT / "studies" / "p4" / "p4b_maxted.json").read_text())
    pts = list(csv.DictReader(open(ROOT / "studies" / "rev2" / "code_sing_points.csv")))
    fig, ax = plt.subplots(1, 3, figsize=(PAGE, 2.6),
                           gridspec_kw=dict(width_ratios=[1, 1, 1.15], wspace=0.42))
    cases = [("phoenix", "none", "PHOENIX"), ("mps2", "none", "MPS-ATLAS 2"),
             ("mps2", "ssd", "+ small-scale\ndynamo"), ("mps2", "B100G", "+ 100 G"),
             ("mps2", "B200G", "+ 200 G")]
    bandmk = (("Kepler", "o"), ("TESS", "s"))              # marker shape = passband
    for a, key, lab, scat in ((ax[0], "dh1", r"$\Delta h'_1$", 0.006),
                              (ax[1], "dh2", r"$\Delta h'_2$", 0.013)):
        a.axvspan(-scat, scat, color=WASH, lw=0, zorder=0)
        a.axvline(0, color=INK2, lw=0.6, zorder=1)
        for k, (band, mk) in enumerate(bandmk):
            for i, (g, s_, _) in enumerate(cases):
                r = [x for x in mx if x["band"] == band and x["grid"] == g and x["state"] == s_
                     and x["law"] == "hybrid5"][0]
                col = BLUE if (g, s_) == ("mps2", "B100G") else INK2
                a.errorbar(r[key][0], i + (k - 0.5) * 0.3, xerr=r[key][1], fmt=mk, ms=3.4,
                           color=col, mfc="white" if band == "TESS" else col, mec=col, mew=0.9,
                           elinewidth=1.0, capsize=0, zorder=3)
        a.set_yticks(range(len(cases)))
        a.set_yticklabels([c[2] for c in cases] if a is ax[0] else [])
        a.set_ylim(len(cases) - 0.5, -0.5)
        a.set_xlabel(lab + "  (observed − predicted)")
        clean_x_grid(a)
    ax[0].set_title("Kepler & TESS hosts · Maxted (2023)")
    jw = [("phoenix|recipe", "PHOENIX", MUTED), ("mps2|recipe", "MPS-ATLAS 2", INK2),
          ("mps2|recipe+B100G", "MPS-ATLAS 2 + 100 G (recommended)", BLUE)]
    stars = [r["star"] for r in pts if r["case"] == "mps2|recipe"]
    ax[2].axvspan(-0.031, 0.031, color=WASH, lw=0, zorder=0)
    ax[2].axvline(0, color=INK2, lw=0.6, zorder=1)
    for k, (case, name, col) in enumerate(jw):
        per = {r["star"]: (float(r["dl"]), float(r["dl_err"])) for r in pts if r["case"] == case}
        v = np.array([per[s_][0] for s_ in stars])
        ax[2].errorbar(v, np.arange(len(stars)) + (k - 1) * 0.24,
                       xerr=[per[s_][1] for s_ in stars], fmt="D", ms=3.2, color=col,
                       mec="white", mew=0.4, elinewidth=1.0, capsize=0, zorder=3)
        ax[2].errorbar(v.mean(), len(stars) + 0.3 + (k - 1) * 0.24,
                       xerr=v.std(ddof=1) / np.sqrt(len(v)), fmt="D", ms=3.6, color=col,
                       mec="white", mew=0.4, elinewidth=1.2, capsize=0, zorder=3)
    ax[2].axhline(len(stars) - 0.2, color=GRID, lw=0.6)
    ax[2].set_yticks(list(range(len(stars))) + [len(stars) + 0.3])
    ax[2].set_yticklabels([("WASP-94A" if s_ == "WASP-94" else s_) for s_ in stars] + ["mean ± SEM"])
    ax[2].set_ylim(len(stars) + 0.8, -0.5)
    ax[2].set_xlabel(r"$\Delta\ell$  (observed − model)")
    ax[2].set_title("JWST hosts · Sing et al. (2026)")
    clean_x_grid(ax[2])
    h1 = [Line2D([], [], ls="", marker=m, ms=4, color=INK2, mfc="white" if b == "TESS" else INK2,
                 mew=0.9, label=b) for b, m in bandmk]
    h2 = [Line2D([], [], ls="", marker="D", ms=4, color=c, label=n) for _, n, c in jw]
    patch = mpl.patches.Patch(color=WASH, label="star-to-star scatter")
    top_legend(fig, h1 + h2 + [patch], y=0.95)
    save(fig, "fig_validation")


def fig_lightcurves():
    """Residual light curves (law − exact model light curve), central transit."""
    d = json.loads((ROOT / "studies" / "lc" / "lc_section.json").read_text())
    p = d["p"]
    t = np.array(d["x"]) / (1 + p)                         # time from mid-transit / (T14/2)
    stars, bnames = ["mid M", "Sun analogue"], d["bands"]
    names = {"mid M": "mid M dwarf (3300 K)", "Sun analogue": "solar analogue (5800 K)"}
    fig, ax = plt.subplots(len(stars), len(bnames), figsize=(PAGE, 3.4), sharex=True, sharey=True,
                           gridspec_kw=dict(hspace=0.34, wspace=0.06))
    for i, st in enumerate(stars):
        for j, bn in enumerate(bnames):
            a = ax[i, j]
            for c in ((1 - p) / (1 + p), 1.0):              # second/third and first/fourth contacts
                for sgn in (-1, 1):
                    a.axvline(sgn * c, color=GRID, lw=0.6, zorder=0)
            a.axhline(0, color=AXIS, lw=0.6, zorder=1)
            for n in ["quadratic", "power2", "hybrid2", "claret4", "hybrid4", "hybrid5"]:
                col, ls, lw, _ = STYLE[n]
                hyb = n.startswith("hybrid")
                a.plot(t, d["resid"][f"{st}|{bn}|{n}"], color=col, ls=ls,
                       lw=1.05 if hyb else 0.75, zorder=4 if n == "hybrid5" else (3 if hyb else 2))
            a.set_yscale("symlog", linthresh=1, linscale=0.5)
            a.set_ylim(-300, 300)
            a.set_yticks([-100, -10, -1, 0, 1, 10, 100])
            a.set_yticklabels(["−100", "−10", "−1", "0", "1", "10", "100"])
            a.set_xlim(-1.08, 1.08)
            w5 = d["worst"][f"{st}|{bn}|hybrid5"]
            wq = d["worst"][f"{st}|{bn}|quadratic"]
            a.text(1.0, 1.03, f"max |ΔF|: hybrid5 {w5:.1f}, quadratic {int(wq + 0.5)} ppm",
                   transform=a.transAxes, fontsize=5.6, color=MUTED, va="bottom", ha="right")
            if i == 0:
                a.text(0.0, 1.22, bn, transform=a.transAxes, fontsize=7.5, fontweight="bold",
                       color=INK, va="bottom")
            if j == 0:
                a.set_ylabel(names[st], fontsize=7.2, fontweight="bold", color=INK)
            if i == len(stars) - 1:
                a.set_xlabel(r"time from mid-transit  /  ($T_{14}/2$)")
    fig.subplots_adjust(left=0.1, top=0.82)
    top_legend(fig, [law_handle(n) for n in LAW_ORDER], y=0.9)
    fig.text(0.005, 0.5, r"law − model  (ppm)", rotation=90, va="center", ha="left",
             fontsize=7.5, color=INK2)
    save(fig, "fig_lightcurves")


if __name__ == "__main__":
    fig_basis()
    fig_lightcurves()
    fig_profiles()
    fig_quality()
    fig_coefficients()
    fig_validation()
