"""One table row: a model atmosphere × passband × law × method, with
coefficients and quality metrics.

    rows = fit_model(spec, "TESS", laws=["hybrid4", "hybrid5"], limb="M1",
                     measure="r", constraint="simplex")
"""

from __future__ import annotations

import numpy as np

from . import bands, fit
from . import laws as L
from . import lightcurve as lc
from .profile import band_profile

CONTROL = ("quadratic", "claret4", "power2")
HYBRIDS = ("hybrid2", "hybrid4", "hybrid5")


def _row(spec, band, profile, name, f, law_eval, lc_metrics):
    m = spec.meta
    row = dict(grid=m.get("grid"), teff=m.get("teff"), logg=m.get("logg"), mh=m.get("mh"),
               band=band, limb=profile.meta.get("limb"), r_c=profile.meta.get("r_c"),
               mu_min_tab=float(profile.mu[0]), n_mu=len(profile.mu),
               law=name, measure=f.measure, constraint=f.constraint, fcm=f.fcm,
               w=[float(x) for x in f.w], rms_I=f.rms, max_abs_I=f.max_abs,
               flux_err=f.flux_err, active=f.active,
               physical=f.extra.get("physical"))
    if lc_metrics:
        e = lc.lc_errors(law_eval, profile)
        row.update(lc_max_ppm=e["max_ppm"], lc_mean_rms_ppm=e["mean_rms_ppm"],
                   lc_worst_by_p=[float(x) for x in e["worst"].max(1)])
    return row


def fit_model(spec, band: str, laws=HYBRIDS + CONTROL, limb: str = "M1",
              measure: str = "r", resample: int | str | None = "exact",
              constraint: str = "simplex", fcm: bool = False, mu_min: float = 0.0,
              lc_metrics: bool = True, r_c: float | None = None) -> list[dict]:
    lam, S = bands.response(band)
    profile = band_profile(spec, lam, S, limb=limb, r_c=r_c)
    rows = []
    for name in laws:
        if name == "power2":
            f = fit.fit_power2(profile, measure, resample, mu_min=mu_min)
            c, a = f.w
            ev = lambda mu, c=c, a=a: L.power2(c, a, mu)
        else:
            law = L.get(name)
            # the simplex is defined for the hybrids only; controls use the exact region
            con = constraint
            if constraint == "simplex" and name not in ("hybrid4", "hybrid5"):
                con = "physical"
            f = fit.fit_linear(law, profile, measure, resample, con, fcm, mu_min)
            ev = lambda mu, law=law, w=f.w: law.intensity(w, mu)
        rows.append(_row(spec, band, profile, name, f, ev, lc_metrics))
    return rows
