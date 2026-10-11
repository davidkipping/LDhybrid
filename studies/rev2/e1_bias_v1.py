"""E1: parameter bias from the limb-darkening law (noiseless).

Truth: the light curve of the PHOENIX band profile (Table 3 cases), planet
p = 0.1 at b = 0.3 on a straight chord z(t) = sqrt(b² + (v t)²), v = 1,
600 points over |t| ≤ 1.2(1 + p). Each law is fitted for (p, b, v):
  fixed — coefficients held at their tabulated (r-weighted) values;
  free  — coefficients fitted too: quadratic and four-term unconstrained;
          hybrid2 in its exact triangle and hybrid4/5 in the simplex, through
          the box coordinates q ∈ [0, 1]^n; power-2 with 0 ≤ c ≤ 1, α > 0.
Reported: Δ(p²) (ppm), Δb, Δv/v, rms residual (ppm).

    .venv/bin/python studies/rev/e1_bias.py        # writes e1_bias.json
"""
import json

import numpy as np
from scipy.optimize import least_squares
from squishierplanet import laws as sl

import rc  # noqa
from ldhybrid import laws, qspace
from ldhybrid import lightcurve as lc

P0, B0, V0 = 0.1, 0.3, 1.0
T = np.linspace(-1.2 * (1 + P0), 1.2 * (1 + P0), 600)


def z_of(b, v):
    return np.sqrt(b * b + (v * T) ** 2)


def blocked_multi(F_of_mu, z, p, n_sub=4):
    """lightcurve.blocked for several profiles at once: F_of_mu(μ) -> (..., k);
    returns ΔF, shape (nz, k). Same quadrature as lightcurve.blocked."""
    z = np.atleast_1d(np.asarray(z, float))
    lo = np.clip(z - p, 0.0, None)
    hi = np.minimum(1.0, z + p)
    bp = np.stack([lo, np.clip(p - z, 0.0, 1.0), np.clip(z - p, 0.0, 1.0),
                   np.minimum(1.0, z + p), hi], axis=1)
    bp = np.sort(np.clip(bp, lo[:, None], hi[:, None]), axis=1)
    edges = [bp[:, :1]]
    for j in range(bp.shape[1] - 1):
        a, b = bp[:, j:j + 1], bp[:, j + 1:j + 2]
        edges.append(a + (b - a) * np.arange(1, n_sub + 1) / n_sub)
    E = np.concatenate(edges, axis=1)
    a, b = E[:, :-1, None], E[:, 1:, None]
    r = a + (b - a) * 0.5 * (1.0 - np.cos(lc._TH))
    dr = (b - a) * 0.5 * np.sin(lc._TH) * lc._WT
    zz = z[:, None, None]
    with np.errstate(divide="ignore", invalid="ignore"):
        cosk = (r * r + zz * zz - p * p) / (2.0 * r * zz)
    kap = np.where(zz == 0, np.where(r < p, np.pi, 0.0),
                   np.arccos(np.clip(np.nan_to_num(cosk, nan=-1.0), -1.0, 1.0)))
    vals = F_of_mu(np.sqrt(np.clip(1.0 - r * r, 0.0, None)))
    return np.einsum("abc,abck->ak", 2.0 * r * kap * dr, vals)


class Model:
    """F(t) for a law at geometry (p, b, v) and coefficients; shape columns
    are cached per geometry so coefficient steps are cheap."""

    def __init__(self, name):
        self.name = name
        self.L = None if name == "power2" else laws.get(name)
        self.cache = {}
        self.fstar = {}

    def _shapes(self, alpha=None):
        if self.name == "power2":
            return lambda mu: np.stack([np.ones_like(mu), 1.0 - mu ** alpha], axis=-1)
        return lambda mu: np.concatenate([np.ones_like(mu)[..., None], self.L.T(mu)], axis=-1)

    def _cols(self, p, b, v, alpha=None):
        k = (p, b, v, alpha)
        if k not in self.cache:
            if len(self.cache) > 64:
                self.cache.clear()
            F = self._shapes(alpha)
            if alpha not in self.fstar:
                th, wt = lc._gl(2000)
                mu = 0.5 * (1.0 - np.cos(th))
                self.fstar[alpha] = 2 * np.pi * np.einsum("i,ik->k", wt * 0.5 * np.sin(th) * mu, F(mu))
            dF = blocked_multi(F, z_of(b, v), p)
            self.cache[k] = [(dF[:, j], self.fstar[alpha][j]) for j in range(dF.shape[1])]
        return self.cache[k]

    def flux(self, p, b, v, w):
        if self.name == "power2":
            c, a = w
            cols = self._cols(p, b, v, a)
            dF = cols[0][0] - c * cols[1][0]
            Fs = cols[0][1] - c * cols[1][1]
        else:
            cols = self._cols(p, b, v)
            dF = cols[0][0] - sum(wj * cj[0] for wj, cj in zip(w, cols[1:]))
            Fs = cols[0][1] - sum(wj * cj[1] for wj, cj in zip(w, cols[1:]))
        return 1.0 - dF / Fs


def coef_map(name, w0):
    """(x0, lower, upper, x -> w) for the free coefficients."""
    n = len(w0)
    if name in ("quadratic", "claret4"):
        return np.array(w0, float), -np.inf * np.ones(n), np.inf * np.ones(n), lambda x: x
    if name == "power2":
        c, a = w0
        return np.array([min(max(c, 1e-6), 1 - 1e-6), a]), np.array([0.0, 1e-4]), \
            np.array([1.0, np.inf]), lambda x: x
    q = np.clip(qspace.to_q(name, w0), 1e-6, 1 - 1e-6)
    if name == "hybrid2":
        f = lambda x: np.asarray(sl.hybrid2_from_q(x[0], x[1], 0.36), float).ravel()
    else:
        f = lambda x: np.asarray(sl.simplex_from_q(x), float)
    return q, np.zeros(n), np.ones(n), f


def run(name, pr, w_tab, F_true, free):
    M = Model(name)
    g0 = np.array([P0, B0, V0])
    if not free:
        res = least_squares(lambda x: (M.flux(*x, w_tab) - F_true) * 1e6, g0,
                            bounds=([0.01, 0.0, 0.2], [0.5, 1.5, 5.0]), x_scale=[0.01, 0.05, 0.05],
                            diff_step=1e-7, xtol=1e-12, ftol=1e-12, gtol=1e-12)
        x, w = res.x, np.array(w_tab, float)
    else:
        c0, lo, hi, fmap = coef_map(name, w_tab)
        res = least_squares(lambda x: (M.flux(*x[:3], fmap(x[3:])) - F_true) * 1e6,
                            np.r_[g0, c0], bounds=(np.r_[[0.01, 0.0, 0.2], lo], np.r_[[0.5, 1.5, 5.0], hi]),
                            diff_step=1e-7, xtol=1e-12, ftol=1e-12, gtol=1e-12)
        x, w = res.x[:3], fmap(res.x[3:])
    p, b, v = x
    rms = float(np.sqrt(np.mean(res.fun ** 2)))
    return dict(dp2_ppm=float((p * p - P0 * P0) * 1e6), db=float(b - B0), dv=float(v / V0 - 1),
                rms_ppm=rms, maxres_ppm=float(np.max(np.abs(res.fun))), w=[float(t) for t in w],
                nfev=int(res.nfev))


def task(star):
    sname, (t, g) = star
    spec = rc.load_t3(t, g)
    out = {}
    for bkey, blab in rc.T3_BANDS:
        pr = rc.profile(spec, bkey)
        F_true = lc.light_curve(pr, z_of(B0, V0), P0)
        for name in rc.LAWS6:
            _, w_tab, _ = rc.fit_law(name, pr, "r")
            out[f"{sname}|{blab}|{name}"] = dict(fixed=run(name, pr, w_tab, F_true, False),
                                                 free=run(name, pr, w_tab, F_true, True))
        print(sname, blab, flush=True)
    return out


if __name__ == "__main__":
    res = rc.run_parallel(task, rc.T3_STARS, workers=5)
    (rc.HERE / "e1_bias.json").write_text(json.dumps({k: v for r in res for k, v in r.items()}))
    print("done")
