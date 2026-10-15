# Inverse q maps: points outside the region must not land in [0, 1]

*For SquishierPlanet, Elómë and Turin. LDhybrid fixed this on 2026-10-15
(commit 272db58).*

## The issue

Every forward map in the family, from q in [0, 1]ⁿ to weights w, has a
square root or another fractional power of q₁ for the depth s = Σw = 1 − I(0):

| Law | Depth |
|---|---|
| Kipping quadratic | s = √q₁ |
| hybrid2 | s = √q₁ |
| hybrid4 | s = q₁^¼ |

Kipping's quadratic map uses s = u₁ + u₂. The inverse computes s from w and
returns q₁ = s² (or s⁴). That power loses the sign of s.

A weight vector with s < 0 has a limb brighter than the centre, so it is
unphysical. Its mirror image −w is often inside the region, and for those
points the inverse returns a q inside [0, 1]ⁿ. The other coordinates are
ratios like t/s, which are unchanged under w → −w, so they stay in range too.

Any code that uses "q in [0, 1]ⁿ" as a physicality check, for example on an
unconstrained fit, a summed coefficient spectrum or a proposal, accepts these
points silently. The same happens with an even power of a negative stick
fraction: in stick-breaking, f^(n−j) with f < 0 and n − j even.

The forward maps are correct, and so are tabulated or sampled values from them.
Only inverse calls on arbitrary w are affected.

## The check

Take the mirror image of points inside the region and invert them. The result
must fall outside the cube:

```python
q = rng.uniform(size=(200, n))
w = -forward(q)                                  # s < 0: unphysical
Q = inverse(w)
assert not np.any(np.all((Q >= 0) & (Q <= 1), axis=-1))
```

In LDhybrid, before the fix, 4,458 of 53,322 unphysical hybrid2 weights from
a uniform box came back with q in [0, 1]², and all of them had s < 0.

## The fix

Keep the sign of every power that can see a negative base. Write sign(x)·|x|ᵖ
instead of xᵖ for even p: q₁ = s·|s| for hybrid2 and the quadratic law, and
sign(s)·s⁴ for hybrid4. The same applies to q₃ = r₂² in hybrid4 and to
f^(n−j) in stick-breaking. In-range points are unchanged.

The alternative is to raise an error outside the region. Elómë's
`hybrid4_to_q_np` already does this, and it is also fine.

## Where it occurs (read-only check of the local checkouts, 2026-10-15)

| Package | Function | Status |
|---|---|---|
| Elómë | `ld.py:u_to_q_np` (quadratic) | `s * s`: mirrors s < 0 |
| Elómë | `ld.py:hybrid2_to_q_np` | `s * s`: mirrors s < 0 |
| Elómë | `ld.py:q_from_simplex_np` (hybrid5) | `f ** (n - j)` with f < 0 is positive for even n − j |
| Elómë | `ld.py:hybrid4_to_q_np` | raises outside the curvature simplex: fine |
| SquishierPlanet | `laws.py:hybrid4_to_q` | `s ** 4` and `(...) ** 2` for q₃ lose the sign |
| SquishierPlanet | hybrid2, simplex | no inverse; forward maps fine |
| Turin | `ldmarg.py:omega_to_q` (quadratic collapsed mode) | `s ** 2`: mirrors s < 0 |
| Turin | `ldhybrid.py:x_to_q` | uses Elómë's raising inverse, then clips to [0, 1]: fine for points on the simplex, but the clip would hide an outside point if the raise were removed |

These are host-side helpers, used to build initial states and convert
outputs, so no posterior is affected while the samplers draw q directly.
They matter when a helper is used to test whether arbitrary weights are
physical.
