"""Figures for the LDhybrid v1.0 paper (static, print).

    .venv/bin/python paper/make_figures.py

Palette: the dataviz reference categorical slots in fixed order (light mode;
validated adjacent-pair for lines, first three all-pairs for scatter/dot).
Laws keep one colour throughout; markers/linestyles double the encoding.
"""

import json
from pathlib import Path

import h5py
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

from ldhybrid import bands, fit, laws
from ldhybrid.grids import phoenix
from ldhybrid.paths import DATA
from ldhybrid.products import constraint_for
from ldhybrid.profile import band_profile

HERE = Path(__file__).parent
FIG = HERE / "figures"
FIG.mkdir(exist_ok=True)
ROOT = HERE.parent

INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
SLOT = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948"]
LAW_ORDER = ["quadratic", "power2", "hybrid2", "claret4", "hybrid4", "hybrid5"]
LAW_COL = dict(zip(LAW_ORDER, SLOT))
LAW_MK = dict(zip(LAW_ORDER, ["o", "s", "^", "D", "v", "P"]))
LAW_LABEL = {"quadratic": "quadratic", "power2": "power-2", "hybrid2": "hybrid2",
             "claret4": "4-term", "hybrid4": "hybrid4", "hybrid5": "hybrid5"}

mpl.rcParams.update({
    "font.size": 9, "axes.labelsize": 9, "axes.titlesize": 9, "legend.fontsize": 8,
    "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
    "text.color": INK, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
    "axes.spines.top": False, "axes.spines.right": False, "lines.linewidth": 1.5,
    "lines.markersize": 4, "legend.frameon": False, "savefig.dpi": 300,
    "savefig.bbox": "tight", "font.family": "serif",
})


def save(fig, name):
    for ext in ("pdf", "png"):
        fig.savefig(FIG / f"{name}.{ext}")
    plt.close(fig)
    print("wrote", name)


def fig_profiles():
    """Model profiles and residuals of the fitted laws, two stars/bands."""
    cases = [((5800, 4.5, "-0.0"), "TESS", "5800 K, TESS"),
             ((3000, 5.0, "-0.0"), "JWST_NIRSpec_G395H", "3000 K, NIRSpec G395H")]
    show = ["quadratic", "claret4", "hybrid4", "hybrid5"]
    fig, ax = plt.subplots(2, 2, figsize=(7.0, 4.2), sharex=True,
                           gridspec_kw=dict(height_ratios=[2, 1.3], hspace=0.08))
    for c, (key, band, title) in enumerate(cases):
        pr = band_profile(phoenix.load(*key), *bands.response(band), limb="M1")
        mu = np.linspace(0, 1, 600)
        y = pr.interpolant()(mu)
        ax[0, c].plot(pr.mu, pr.y, "o", color=INK2, ms=3, label="PHOENIX nodes", zorder=3)
        ax[0, c].plot(mu, y, color=INK2, lw=0.8)
        for n in show:
            L = laws.get(n)
            w = fit.fit_linear(L, pr, "r", 2000, constraint_for(n)).w
            ax[1, c].plot(mu, (L.intensity(w, mu) - y) * 1e3, color=LAW_COL[n], label=LAW_LABEL[n])
        ax[0, c].set_title(title, loc="left")
        ax[1, c].axhline(0, color=INK2, lw=0.6)
        ax[1, c].set_xlabel(r"$\mu$")
        ax[1, c].set_yscale("symlog", linthresh=2)
        ax[1, c].set_ylim(-300, 300)
    ax[0, 0].set_ylabel(r"$I(\mu)/I(1)$")
    ax[1, 0].set_ylabel(r"law $-$ model ($10^{-3}$)")
    h, l = ax[1, 0].get_legend_handles_labels()
    fig.legend(h, l, ncol=4, loc="upper center", bbox_to_anchor=(0.5, 0.0))
    ax[0, 0].legend(loc="lower right")
    save(fig, "fig_profiles")


def fig_quality():
    """Median and 99th-percentile light-curve error vs Teff, all bands, per law."""
    fig, ax = plt.subplots(1, 2, figsize=(7.0, 2.9), sharey=True)
    for a, grid, title in ((ax[0], "phoenix", "PHOENIX v3.0"), (ax[1], "mps2", "MPS-ATLAS Set 2")):
        h = h5py.File(DATA / "products" / f"{grid}.h5")
        T, cov = h["teff"][()], h["covered"][()]
        Tb = np.unique(np.round(T / 200) * 200)
        for n in LAW_ORDER:
            e = h[f"{n}|r|lc_ppm"][()].astype(float)
            e[~cov] = np.nan
            med = [np.nanmedian(e[np.abs(T - t) <= 100]) for t in Tb]
            a.plot(Tb, med, color=LAW_COL[n], marker=LAW_MK[n], ms=3.5, label=LAW_LABEL[n])
        a.set_yscale("log")
        a.set_title(title, loc="left")
        a.set_xlabel(r"$T_{\rm eff}$ (K)")
    ax[0].set_ylabel("median light-curve error (ppm)")
    h, l = ax[0].get_legend_handles_labels()
    fig.legend(h, l, ncol=6, loc="upper center", bbox_to_anchor=(0.5, -0.04))
    save(fig, "fig_quality")


def fig_coefficients():
    """hybrid5 weights vs Teff (log g 4.5, [M/H] = 0) for three bands, PHOENIX
    (solid) and MPS-ATLAS Set 2 (dashed)."""
    bset = [("TESS", "TESS"), ("HST_WFC3_G141", "WFC3 G141"), ("JWST_NIRSpec_G395H", "NIRSpec G395H")]
    fig, ax = plt.subplots(1, 5, figsize=(7.2, 2.2), sharex=True)
    for grid, ls in (("phoenix", "-"), ("mps2", "--")):
        h = h5py.File(DATA / "products" / f"{grid}.h5")
        B = [b.decode() for b in h["bands"][()]]
        sel = (h["logg"][()] == 4.5) & (np.abs(h["mh"][()]) < 1e-9)
        T = h["teff"][()][sel]
        o = np.argsort(T)
        W = h["hybrid5|r|w"][()][sel][o]
        for k, (b, lab) in enumerate(bset):
            for j in range(5):
                ax[j].plot(T[o], W[:, B.index(b), j], ls=ls, color=SLOT[k],
                           label=lab if (grid == "phoenix" and j == 0) else None)
    for j in range(5):
        ax[j].set_title(f"$w_{j + 1}$", loc="left")
        ax[j].set_xlabel(r"$T_{\rm eff}$ (K)")
        ax[j].tick_params(axis="x", labelsize=7)
    from matplotlib.lines import Line2D
    h, l = ax[0].get_legend_handles_labels()
    h += [Line2D([], [], color=INK2, ls="-"), Line2D([], [], color=INK2, ls="--")]
    l += ["PHOENIX v3.0", "MPS-ATLAS Set 2"]
    fig.tight_layout()
    fig.legend(h, l, ncol=5, loc="upper center", bbox_to_anchor=(0.5, 0.0))
    save(fig, "fig_coefficients")


def fig_validation():
    """Kepler h'1 offsets (Maxted 2023) and JWST ℓ offsets (Sing et al. 2026)."""
    mx = json.loads((ROOT / "studies" / "p4" / "p4b_maxted.json").read_text())
    sg = json.loads((ROOT / "studies" / "p4" / "p4d_sing.json").read_text())
    fig, ax = plt.subplots(1, 3, figsize=(7.2, 2.6), gridspec_kw=dict(width_ratios=[1, 1, 1.2]))
    cases = [("phoenix", "none", "PHOENIX"), ("mps2", "none", "MPS-ATLAS 2"),
             ("mps2", "ssd", "MPS2 + ssd"), ("mps2", "B100G", "MPS2 + 100 G"),
             ("mps2", "B200G", "MPS2 + 200 G")]
    for a, key, lab in ((ax[0], "dh1", r"$\Delta h'_1$"), (ax[1], "dh2", r"$\Delta h'_2$")):
        for k, (band, col) in enumerate((("Kepler", SLOT[0]), ("TESS", SLOT[1]))):
            for i, (g, s, name) in enumerate(cases):
                r = [x for x in mx if x["band"] == band and x["grid"] == g and x["state"] == s
                     and x["law"] == "hybrid5"][0]
                a.errorbar(r[key][0], i + (k - 0.5) * 0.25, xerr=r[key][1], fmt="o", color=col,
                           ms=4, capsize=0, label=band if i == 0 else None)
        a.axvline(0, color=INK2, lw=0.6)
        a.set_yticks(range(len(cases)))
        a.set_yticklabels([c[2] for c in cases] if a is ax[0] else [])
        a.invert_yaxis()
        a.set_xlabel(lab + " (obs $-$ pred)")
    ax[0].legend(loc="upper center", bbox_to_anchor=(1.0, -0.28), ncol=2, fontsize=7)
    ax[0].set_title("Kepler/TESS, Maxted (2023)", loc="left")
    jw = [("phoenix|hybrid5", "PHOENIX", SLOT[0]), ("mps2|hybrid5", "MPS-ATLAS 2", SLOT[1]),
          ("mps2|hybrid5+B100G", "MPS2 + 100 G", SLOT[2])]
    stars = [s for s, _, _ in sg["mps2|hybrid5"]["per_star"]]
    for k, (case, name, col) in enumerate(jw):
        per = {s: v for s, v, _ in sg[case]["per_star"]}
        ax[2].plot([per[s] for s in stars], np.arange(len(stars)) + (k - 1) * 0.22, "o",
                   color=col, ms=4, label=name)
    ax[2].axvline(0, color=INK2, lw=0.6)
    ax[2].set_yticks(range(len(stars)))
    ax[2].set_yticklabels(stars, fontsize=7)
    ax[2].invert_yaxis()
    ax[2].set_xlabel(r"$\Delta\ell$ (obs $-$ model)")
    ax[2].set_title("JWST, Sing et al. (2026)", loc="left")
    ax[2].legend(loc="upper center", bbox_to_anchor=(0.45, -0.28), ncol=3, fontsize=7)
    fig.tight_layout()
    save(fig, "fig_validation")


if __name__ == "__main__":
    fig_profiles()
    fig_quality()
    fig_coefficients()
    fig_validation()
