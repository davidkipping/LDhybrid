"""Gate G1(b): reproduce ExoTiC-LD (v3.2) coefficients on the same intensity
files (its Kurucz and Stagger bundle), with our band integration and fitting.

ExoTiC-LD's recipe: photon-unit intensities and throughputs, linearly
interpolated and integrated by 500-node Gauss–Legendre over the mode's range;
I(μ=1) normalisation; unweighted Levenberg–Marquardt on the native μ nodes
with μ ≥ 0.1. Ours: exact bin integrals, then measure="nodes", mu_min=0.1.

    .venv/bin/python studies/g1/g1_exoticld.py
"""

import json
import os
from pathlib import Path

import numpy as np
from exotic_ld import StellarLimbDarkening

from ldhybrid import bands, fit, laws
from ldhybrid.grids import exoticld
from ldhybrid.profile import Profile

HERE = Path(__file__).parent
LINK = exoticld.ROOT / "Sensitivity_files"
MODELS = {"kurucz": [(0.0, 4500, 4.5), (0.0, 5750, 4.5), (-1.0, 6500, 4.0), (0.0, 3500, 5.0)],
          "stagger": [(0.0, 4500, 4.5), (0.0, 5500, 4.5), (-1.0, 6500, 4.0), (0.0, 6000, 4.5)]}
MODES = {"TESS": "TESS", "HST_WFC3_G141": "HST_WFC3_G141", "HST_STIS_G430L": "HST_STIS_G430L",
         "JWST_NIRSpec_PRISM": "JWST_NIRSpec_prism", "JWST_NIRISS_SOSS_o1": "JWST_NIRISS_SOSSo1",
         "JWST_NIRSpec_G395H": "JWST_NIRSpec_G395H", "Spitzer_IRAC2": "Spitzer_IRAC_Ch2"}


def main():
    if not LINK.exists():
        os.symlink(os.path.relpath(exoticld.RAW / "bands" / "exoticld", exoticld.ROOT), LINK)
    rows = []
    for grid, models in MODELS.items():
        for mh, teff, logg in models:
            spec = exoticld.load(grid, mh, teff, logg)
            sld = StellarLimbDarkening(mh, teff, logg, grid, str(exoticld.ROOT), verbose=0)
            for ours, theirs in MODES.items():
                lam, S = bands.response(ours)
                wr = [lam[S > 0].min(), lam[S > 0].max()]
                q_t = np.array(sld.compute_quadratic_ld_coeffs(wr, theirs, mu_min=0.1))
                c_t = np.array(sld.compute_4_parameter_non_linear_ld_coeffs(wr, theirs, mu_min=0.1))
                I_t = sld.I_mu.copy()                              # their band profile
                I_o = spec.band(lam, S, min_coverage=0.99)
                pr = Profile(spec.mu, I_o)
                q_o = fit.fit_linear(laws.quadratic(), pr, "nodes", mu_min=0.1).w
                c_o = fit.fit_linear(laws.claret4(), pr, "nodes", mu_min=0.1).w
                # same fit on *their* profile: isolates fitting from integration
                pr_t = Profile(sld.mus, I_t)
                c_ot = fit.fit_linear(laws.claret4(), pr_t, "nodes", mu_min=0.1).w
                y_t = np.sort(I_t / I_t[0])
                rows.append(dict(grid=grid, mh=mh, teff=teff, logg=logg, band=ours,
                                 dprofile=float(np.max(np.abs(np.sort(pr.y) - y_t))),
                                 dquad=float(np.max(np.abs(q_o - q_t))),
                                 dclaret4=float(np.max(np.abs(c_o - c_t))),
                                 dclaret4_same_profile=float(np.max(np.abs(c_ot - c_t)))))
                r = rows[-1]
                print(f"{grid:7s} {mh:+.1f} {teff} {logg} {ours:22s} |dI| {r['dprofile']:.1e} "
                      f"|dquad| {r['dquad']:.1e} |d4| {r['dclaret4']:.1e} "
                      f"(same profile {r['dclaret4_same_profile']:.1e})", flush=True)
    (HERE / "g1_exoticld.json").write_text(json.dumps(rows, indent=1))
    for k in ("dprofile", "dquad", "dclaret4", "dclaret4_same_profile"):
        print(f"max {k}: {max(r[k] for r in rows):.2e}")


if __name__ == "__main__":
    main()
