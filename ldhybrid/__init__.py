"""LDhybrid: theoretical limb-darkening coefficient tables for the rational
hybrid2/4/5 laws of SquishierPlanet. See docs/plan.md and docs/products.md."""

__version__ = "1.0.0.dev0"

from .api import (GRIDS, coefficients, coefficients_custom, condition,  # noqa: E402,F401
                  observables, prior, to_law)
