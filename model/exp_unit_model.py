"""Golden model for rtl/exp_unit.sv.

e^x = 2^(x*log2(e)) = 2^expo * 2^frac, with expo = round(t) and frac = t - expo
in [-1/2, 1/2). The exponent costs a shift; only 2^frac needs a polynomial.

This mirrors the RTL stage for stage. It is not a reference exp() -- it is a
model of the circuit, which is what makes bit-exact comparison the right test.
"""

import numpy as np

import poly_eval_model as poly
from fixedpoint import Fmt, Fx

LOG2E = 1.4426950408889634


class ExpFormats:
    """The parameters of exp_unit.sv, and the shifts derived from them."""

    def __init__(self, X_W=18, X_F=17, L_W=18, L_F=16, T_W=20, T_F=18,
                 C_W=18, C_F=16, W_W=18, W_F=16, Y_W=32, Y_F=31):
        self.x = Fmt(X_W, X_F, "half_up")       # module input
        self.l = Fmt(L_W, L_F, "half_up")       # the log2(e) constant
        self.t = Fmt(T_W, T_F, "half_up")       # t + 1/2, the value that gets split
        self.arg = Fmt(T_F, T_F, "half_up")     # frac into poly_eval: width == frac
        self.c = Fmt(C_W, C_F, "half_up")       # coefficients
        self.w = Fmt(W_W, W_F, "half_up")       # poly_eval working format
        self.T_F, self.W_F, self.Y_W, self.Y_F = T_F, W_F, Y_W, Y_F
        self.SH_Y = Y_F - W_F
        self.EXPO_W = T_W - T_F

    @property
    def max_shift(self):
        """Largest right shift the barrel shifter must do, from x's lower bound."""
        x_min = -(2.0 ** (self.x.width - 1 - self.x.frac))
        return int(-np.floor(x_min * LOG2E + 0.5))


DEFAULT = ExpFormats()

# The wide instance, for norm_cdf's density: x down to -16.
DENSITY = ExpFormats(X_W=24, X_F=19, L_W=24, L_F=22, T_W=24, T_F=18,
                     C_W=18, C_F=16, W_W=18, W_F=16, Y_W=32, Y_F=31)


def expo_boundaries(fmt=DEFAULT, span=3):
    """x values straddling every expo transition, span LSBs either side.

    expo = floor(t + 1/2) changes wherever t = k - 1/2, so x = (k - 1/2)/log2(e).
    Across that seam expo jumps by one while frac snaps from +1/2 to -1/2, and
    2^expo * 2^frac has to come out continuous. A wrong bias, a missed MSB flip
    or a bad sign extension on expo is correct everywhere except within an LSB
    of these points, which a uniform sweep steps straight over.
    """
    lsb = 2.0 ** -fmt.x.frac
    x_min = -(2.0 ** (fmt.x.width - 1 - fmt.x.frac))
    edges = (np.arange(-fmt.max_shift, 1) - 0.5) / LOG2E
    offsets = np.arange(-span, span + 1) * lsb
    xs = (edges[:, None] + offsets[None, :]).ravel()
    return xs[(xs >= x_min) & (xs <= 0.0)]


def stages(x, fmt=DEFAULT, coeffs=None):
    """Every intermediate, keyed by the signal name in exp_unit.sv."""
    coeffs = poly.exp2_coefficients() if coeffs is None else coeffs

    # clocks 0-1. A product is exact: widths add, nothing rounds.
    x_r = fmt.x.const(x)
    t_prod = x_r * fmt.l.const(LOG2E)

    # clock 2. Narrow to T_F fractional bits, round half up, and add the 1/2
    # bias. In the RTL all three ride on the single constant RND_BIAS.
    t_biased = fmt.t.cast(t_prod + Fx(1, 1))

    # The split is wiring. Python's >> on a negative int is already floor, which
    # is what an arithmetic shift gives and what round(t) = floor(t + 1/2) needs.
    expo = t_biased.code >> fmt.T_F
    low = t_biased.code & ((1 << fmt.T_F) - 1)
    # low holds frac + 1/2 in excess-1/2 form; subtracting half is the MSB flip.
    frac = Fx(low - (1 << (fmt.T_F - 1)), fmt.T_F)

    # clocks 3-9. frac is exactly representable in the argument format, so the
    # float roundtrip through const() cannot re-round it.
    signif = poly.evaluate_estrin_fixed(coeffs, frac.to_float(), fmt.c, fmt.arg, fmt.w)

    # clocks 10-11. y = signif * 2^expo. The left align is constant (wiring) and
    # carries one spare bit; the right shift is the barrel shifter; the final
    # (+1) >> 1 rounds half up at a fixed position and drops the guard.
    rshift = -expo
    signif_algn = signif.code << (fmt.SH_Y + 1)
    sh1 = signif_algn >> rshift
    y = (sh1 + 1) >> 1

    return {
        "x_r": x_r, "t_prod": t_prod, "t_biased": t_biased,
        "expo": expo, "frac": frac, "signif": signif,
        "signif_algn": signif_algn, "rshift": rshift, "sh1": sh1,
        "y": y & ((1 << fmt.Y_W) - 1),
    }


def evaluate_fixed(x, fmt=DEFAULT, coeffs=None):
    """The y codes the RTL should produce, for the same x."""
    return stages(x, fmt, coeffs)["y"]


def evaluate_float(x, fmt=DEFAULT, coeffs=None):
    """Same, as a real number."""
    return np.asarray(evaluate_fixed(x, fmt, coeffs), dtype=float) * 2.0 ** -fmt.Y_F


def relative_error(x, fmt=DEFAULT, coeffs=None):
    ref = np.exp(np.asarray(x, dtype=float))
    return (evaluate_float(x, fmt, coeffs) - ref) / ref


def sweep(fmt=DEFAULT, count=20001, coeffs=None):
    """Relative accuracy over the whole contract, overall and by decade of y."""
    x_min = -(2.0 ** (fmt.x.width - 1 - fmt.x.frac))
    xs = np.linspace(x_min, 0.0, count)
    rel = np.abs(relative_error(xs, fmt, coeffs))
    bits = lambda r: float("inf") if r <= 0 else -np.log2(r)

    print(f"x in [{x_min:g}, 0], {count} points, max_shift = {fmt.max_shift}")
    print(f"  overall  max rel {rel.max():.3e}  ({bits(rel.max()):.1f} bits)"
          f"   rms {np.sqrt((rel ** 2).mean()):.3e}")
    edges = np.linspace(x_min, 0.0, 5)
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (xs >= lo) & (xs <= hi)
        print(f"  x in [{lo:7.2f},{hi:7.2f}]  max rel {rel[m].max():.3e}"
              f"  ({bits(rel[m].max()):.1f} bits)")
    return rel


def trace(x, fmt=DEFAULT, coeffs=None):
    """One value, every intermediate printed. For chasing an RTL mismatch."""
    s = stages(np.asarray([x], dtype=float), fmt, coeffs)
    code = lambda v: int(v.code[0]) if isinstance(v, Fx) else int(np.asarray(v)[0])
    real = lambda v: code(v) * 2.0 ** -v.frac if isinstance(v, Fx) else code(v)

    print(f"x           = {x:+.9f}")
    for name in ("x_r", "t_prod", "t_biased", "expo", "frac", "signif",
                 "signif_algn", "rshift", "sh1", "y"):
        v = s[name]
        suffix = f"   = {real(v):+.9f}" if isinstance(v, Fx) else ""
        print(f"  {name:12s} {code(v):>16d}{suffix}")

    got = code(s["y"]) * 2.0 ** -fmt.Y_F
    ref = float(np.exp(x))
    rel = (got - ref) / ref
    print(f"  y           = {got:.12f}   exp(x) = {ref:.12f}")
    print(f"  rel {rel:+.3e}  ({-np.log2(abs(rel)):.1f} bits)")


if __name__ == "__main__":
    trace(-0.75)
    print()
    sweep()
    print()
    sweep(DENSITY)
