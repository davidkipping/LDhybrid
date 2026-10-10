# Draft: request for PHOENIX NewEra centre-to-limb data (not sent)

To: Peter Hauschildt (Hamburger Sternwarte) — address to confirm
Subject: Request for NewEra specific-intensity (centre-to-limb) data

Dear Peter,

I am building theoretical limb-darkening tables for a new family of rational
limb-darkening laws. They give elementary transit light curves for oblate
planets, and the coefficients can be fitted as a convex problem, so positivity
and monotonicity hold by construction. Our baseline is the
PHOENIX-ACES SpecInt v3.0 grid, which stops at 5.5 µm.

Hauschildt et al. (2025, A&A 698, A47) note that NewEra limb-darkening data are
available on request. Would you be willing to share the NewEra specific
intensities I(λ, μ) for FGKM dwarfs and subgiants? Our scope is roughly
Teff 2300–7500 K, log g 3.5–5.5 and [M/H] −2 to +1. We are most interested in:

* the full μ sampling (127 points?), including the outer layers, so that we
  can apply and compare limb definitions (max |dI/dr|, Wittkowski M1/M2);
* wavelength coverage to ≥ 14 µm, for JWST MIRI LRS, especially for M dwarfs;
* any resolution you can share. Photon-binned at R ≈ 3000 would be enough.

We will cite the NewEra papers and acknowledge the data, and we are happy to
share our coefficient tables and code before publication.

With best wishes,
David Kipping
