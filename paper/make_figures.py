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

from ldhybrid import bands, fit, laws
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
    show = ["quadratic", "claret4", "hybrid4", "hybrid5"]
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
            L = laws.get(n)
            w = fit.fit_linear(L, pr, "r", 2000, constraint_for(n)).w
            col, ls, lw, _ = STYLE[n]
            top.plot(mu, L.intensity(w, mu), color=col, ls=ls, lw=lw * 0.8, zorder=3)
            ins.plot(mu, L.intensity(w, mu), color=col, ls=ls, lw=lw, zorder=3)
            bot.plot(mu, (L.intensity(w, mu) - y) * 1e3, color=col, ls=ls, lw=lw)
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
    top_legend(fig, [nodes] + [law_handle(n) for n in show[::-1]], y=0.945)
    save(fig, "fig_profiles")


def fig_quality():
    fig, ax = plt.subplots(1, 2, figsize=(PAGE, 2.6), sharey=True, gridspec_kw=dict(wspace=0.08))
    for a, grid, title in ((ax[0], "phoenix", "PHOENIX v3.0 (spherical)"),
                           (ax[1], "mps2", "MPS-ATLAS Set 2 (plane-parallel)")):
        h = h5py.File(DATA / "products" / f"{grid}.h5")
        T, cov = h["teff"][()], h["covered"][()]
        Tb = np.unique(np.round(T / 200) * 200)
        for n in LAW_ORDER[::-1]:
            e = h[f"{n}|r|lc_ppm"][()].astype(float)
            e[~cov] = np.nan
            sel = [np.abs(T - t) <= 100 for t in Tb]
            med = np.array([np.nanmedian(e[s]) for s in sel])
            col, ls, lw, _ = STYLE[n]
            if n.startswith("hybrid"):
                lo = np.array([np.nanpercentile(e[s], 16) for s in sel])
                hi = np.array([np.nanpercentile(e[s], 84) for s in sel])
                a.fill_between(Tb, lo, hi, color=col, alpha=0.13, lw=0)
            a.plot(Tb, med, color=col, ls=ls, lw=lw)
        a.set_yscale("log")
        a.set_ylim(0.3, 300)
        a.set_title(title)
        a.set_xlabel(r"$T_{\rm eff}$ (K)")
    ax[0].set_ylabel("light-curve error (ppm)")
    ax[0].yaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda v, _: f"{v:g}"))
    top_legend(fig, [law_handle(n) for n in LAW_ORDER], y=0.95)
    fig.text(0.995, -0.02, "lines: median over passbands, gravities and metallicities;"
             " bands: 16th–84th percentile (hybrid laws)", ha="right", va="top",
             fontsize=6, color=MUTED)
    save(fig, "fig_quality")


def fig_coefficients():
    bset = [("Kepler", "Kepler", 250), ("TESS", "TESS", 350),
            ("HST_WFC3_G141", "WFC3 G141", 450), ("JWST_NIRISS_SOSS_o1", "NIRISS SOSS", 550),
            ("JWST_NIRSpec_G395H", "NIRSpec G395H", 650)]
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
                ax[r, j].plot(T[o], W[:, B.index(b), j], color=RAMP[step], lw=1.3)
        ax[r, 0].annotate(label, xy=(-0.42, 0.5), xycoords="axes fraction", rotation=90,
                          ha="center", va="center", fontsize=7.5, fontweight="bold", color=INK)
    for j in range(5):
        ax[0, j].set_title(f"$w_{j + 1}$")
        ax[1, j].set_xlabel(r"$T_{\rm eff}$ (K)")
        ax[1, j].set_xlim(2300, 7500)
        ax[1, j].xaxis.set_major_locator(mpl.ticker.MultipleLocator(2000))
    handles = [Line2D([], [], color=RAMP[s], lw=1.6, label=lab) for _, lab, s in bset]
    top_legend(fig, handles, y=0.955)
    fig.text(0.995, -0.015, r"hybrid5, $\log g = 4.5$, [M/H] = 0; passbands ordered by wavelength",
             ha="right", va="top", fontsize=6, color=MUTED)
    save(fig, "fig_coefficients")


def fig_validation():
    mx = json.loads((ROOT / "studies" / "p4" / "p4b_maxted.json").read_text())
    sg = json.loads((ROOT / "studies" / "p4" / "p4d_sing.json").read_text())
    fig, ax = plt.subplots(1, 3, figsize=(PAGE, 2.5),
                           gridspec_kw=dict(width_ratios=[1, 1, 1.15], wspace=0.42))
    cases = [("phoenix", "none", "PHOENIX"), ("mps2", "none", "MPS-ATLAS 2"),
             ("mps2", "ssd", "+ small-scale dynamo"), ("mps2", "B100G", "+ 100 G"),
             ("mps2", "B200G", "+ 200 G")]
    bandmk = (("Kepler", INK2), ("TESS", "white"))           # filled / open markers
    for a, key, lab, scat in ((ax[0], "dh1", r"$\Delta h'_1$", 0.006),
                              (ax[1], "dh2", r"$\Delta h'_2$", 0.013)):
        a.axvspan(-scat, scat, color=WASH, lw=0, zorder=0)
        a.axvline(0, color=INK2, lw=0.6, zorder=1)
        for k, (band, face) in enumerate(bandmk):
            for i, (g, s, _) in enumerate(cases):
                r = [x for x in mx if x["band"] == band and x["grid"] == g and x["state"] == s
                     and x["law"] == "hybrid5"][0]
                col = BLUE if (g, s) == ("mps2", "B100G") else INK2
                a.errorbar(r[key][0], i + (k - 0.5) * 0.3, xerr=r[key][1], fmt="o", ms=3.6,
                           color=col, mfc=col if face != "white" else "white", mec=col, mew=0.9,
                           elinewidth=1.0, capsize=0, zorder=3)
        a.set_yticks(range(len(cases)))
        a.set_yticklabels([c[2] for c in cases] if a is ax[0] else [])
        a.set_ylim(len(cases) - 0.5, -0.5)
        a.set_xlabel(lab + "  (observed − predicted)")
        clean_x_grid(a)
    ax[0].set_title("Kepler & TESS hosts · Maxted (2023)")
    jw = [("phoenix|hybrid5", "PHOENIX", MUTED), ("mps2|hybrid5", "MPS-ATLAS 2", INK2),
          ("mps2|hybrid5+B100G", "MPS-ATLAS 2 + 100 G (recommended)", BLUE)]
    stars = [s for s, _, _ in sg["mps2|hybrid5"]["per_star"]]
    ax[2].axvspan(-0.03, 0.03, color=WASH, lw=0, zorder=0)
    ax[2].axvline(0, color=INK2, lw=0.6, zorder=1)
    for k, (case, name, col) in enumerate(jw):
        per = {s: v for s, v, _ in sg[case]["per_star"]}
        ax[2].plot([per[s] for s in stars], np.arange(len(stars)) + (k - 1) * 0.24, "o", ms=3.6,
                   color=col, mec="white", mew=0.5, zorder=3)
    ax[2].set_yticks(range(len(stars)))
    ax[2].set_yticklabels(stars)
    ax[2].set_ylim(len(stars) - 0.5, -0.5)
    ax[2].set_xlabel(r"$\Delta\ell$  (observed − model)")
    ax[2].set_title("JWST hosts · Sing et al. (2026)")
    clean_x_grid(ax[2])
    h1 = [Line2D([], [], ls="", marker="o", ms=4, color=INK2, mfc=f if f == "white" else INK2,
                 mew=0.9, label=b) for b, f in bandmk]
    h2 = [Line2D([], [], ls="", marker="o", ms=4, color=c, label=n) for _, n, c in jw]
    patch = mpl.patches.Patch(color=WASH, label="star-to-star scatter")
    top_legend(fig, h1 + h2 + [patch], y=0.95)
    save(fig, "fig_validation")


if __name__ == "__main__":
    fig_profiles()
    fig_quality()
    fig_coefficients()
    fig_validation()
