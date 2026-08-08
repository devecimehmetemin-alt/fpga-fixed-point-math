"""Golden model for rtl/ln_unit.sv.

ln x = e*ln2 + ln(1+z), for x = m * 2^e out of lzc_norm and z = m - 1 in [0,1).
m's leading bit is the integer 1, so dropping it *is* the subtraction: the range
reduction is a slice and nothing else.

x arrives as a raw W-bit code read as Q(W-F).F, so the exponent the caller means
is e - F. y is signed with W_F fractional bits.

This mirrors the RTL stage for stage, not a reference log(). The two places that
matter: z truncates rather than rounds, and the back-end add happens at W_F with
no realignment because LN2 is stored at W_F too.
"""

import numpy as np

import lzc_norm_model as lzc
import poly_eval_model as poly
from fixedpoint import Fmt, Fx

LN2 = 0.6931471805599453


class LnFormats:
    """The parameters of ln_unit.sv, and the widths derived from them."""

    def __init__(self, W=32, F=28, X_W=18, X_F=17, C_W=18, C_F=17,
                 W_W=18, W_F=16, L_W=18, Y_W=24):
        self.c = Fmt(C_W, C_F, "half_up") # coefficients
        self.arg = Fmt(X_W, X_F, "half_up") # z into poly_eval
        self.w = Fmt(W_W, W_F, "half_up") # poly_eval working format
        self.ln2 = Fmt(L_W, W_F, "half_up") # the constant, at the working scale
        self.W, self.F, self.X_F, self.W_F, self.Y_W = W, F, X_F, W_F, Y_W

    @property
    def ln2_code(self):
        return int(self.ln2.const(LN2).code)


DEFAULT = LnFormats()


def _wrap(code, width):
    # Y_W'(...) in the RTL truncates, so wrap the same way and stay bit-exact.
    # y overflowing is a parameter bug rather than a rounding effect, though,
    # so it is worth saying out loud when it happens.
    span = 1 << width
    wrapped = ((code + (span >> 1)) % span) - (span >> 1)
    if np.any(wrapped != code):
        print(f"warning: y does not fit in Y_W={width} bits and wrapped")
    return wrapped


def stages(x, fmt=DEFAULT, coeffs=None):
    """Every intermediate, keyed by the signal name in ln_unit.sv.

    x is raw W-bit unsigned integer codes, not reals.
    """
    coeffs = poly.ln1p_coefficients() if coeffs is None else coeffs
    x = np.atleast_1d(np.asarray(x, dtype=np.int64))

    # clocks 0-1, lzc_norm
    norm = [lzc.normalize(int(v), fmt.W) for v in x]
    m = np.array([n[0] for n in norm], dtype=np.int64)
    e = np.array([n[1] for n in norm], dtype=np.int64)
    zero = np.array([n[3] for n in norm], dtype=bool)

    # clock 2. m[W-2 -: X_F] keeps X_F bits from under the leading one. This
    # truncates, and at X_F = 17 the 7.6e-06 it throws away is the same size as
    # the polynomial's own error -- which is why X_F was picked there.
    z = Fx((m >> (fmt.W - 1 - fmt.X_F)) & ((1 << fmt.X_F) - 1), fmt.X_F)

    # F comes off once, here. lzc_norm reports the exponent of x read as an
    # integer; a Q(W-F).F caller means e - F, and skipping it is a silent 2^F.
    et = e - fmt.F

    # clocks 3-9, ln(1+z)
    l = poly.evaluate_estrin_fixed(coeffs, z.to_float(), fmt.c, fmt.arg, fmt.w)

    # clock 9. et is an integer and LN2 already carries W_F fractional bits, so
    # the product lands at W_F and the add below needs no shift.
    prod = et * fmt.ln2_code

    # clock 10
    y = _wrap(np.where(zero, 0, prod + l.code), fmt.Y_W)

    return {
        "m": m, "e": e, "zero": zero, "z": z, "e_d": et,
        "l": l, "prod": prod, "y": y,
    }


def evaluate_fixed(x, fmt=DEFAULT, coeffs=None):
    """The y codes the RTL should produce, for the same x codes."""
    return stages(x, fmt, coeffs)["y"]


def evaluate_float(x, fmt=DEFAULT, coeffs=None):
    """Same, as a real number."""
    return np.asarray(evaluate_fixed(x, fmt, coeffs), dtype=float) * 2.0 ** -fmt.W_F


def absolute_error(x, fmt=DEFAULT, coeffs=None):
    """y - ln(x / 2^F). Absolute, deliberately: d(ln x) = dx/x, so absolute
    error in a log is relative error in its argument, which is what the caller
    downstream actually cares about. Relative error of ln is unbounded at x = 1
    and says nothing useful."""
    ref = np.log(np.asarray(x, dtype=float) * 2.0 ** -fmt.F)
    return evaluate_float(x, fmt, coeffs) - ref


def sweep(fmt=DEFAULT, count=20001, coeffs=None):
    """Absolute accuracy over the whole domain, overall and by octave of x."""
    hi = (1 << fmt.W) - 1
    xs = np.unique(np.exp(np.linspace(0.0, np.log(hi), count)).round().astype(np.int64))
    xs = xs[xs > 0]
    err = np.abs(absolute_error(xs, fmt, coeffs))
    real = xs * 2.0 ** -fmt.F
    bits = lambda r: float("inf") if r <= 0 else -np.log2(r)

    print(f"x in [{real[0]:.3e}, {real[-1]:g}], {len(xs)} codes")
    print(f"  overall  max abs {err.max():.3e}  ({bits(err.max()):.1f} bits)"
          f"   rms {np.sqrt((err ** 2).mean()):.3e}")
    edges = np.geomspace(real[0], real[-1], 5)
    for lo, up in zip(edges[:-1], edges[1:]):
        sel = (real >= lo) & (real <= up)
        print(f"  x in [{lo:9.3e},{up:9.3e}]  max abs {err[sel].max():.3e}"
              f"  ({bits(err[sel].max()):.1f} bits)")
    return err


def trace(x, fmt=DEFAULT, coeffs=None):
    """One code, every intermediate printed. For chasing an RTL mismatch."""
    s = stages([x], fmt, coeffs)
    code = lambda v: int(v.code[0]) if isinstance(v, Fx) else int(np.atleast_1d(v)[0])

    print(f"x = {x} = {x * 2.0 ** -fmt.F:.9f}")
    for name in ("m", "e", "zero", "z", "e_d", "l", "prod", "y"):
        v = s[name]
        suffix = ""
        if name == "z":
            suffix = f"   = {code(v) * 2.0 ** -fmt.X_F:+.9f}"
        elif name in ("l", "prod", "y"):
            suffix = f"   = {code(v) * 2.0 ** -fmt.W_F:+.9f}"
        print(f"  {name:6s} {code(v):>14d}{suffix}")

    got = code(s["y"]) * 2.0 ** -fmt.W_F
    ref = float(np.log(x * 2.0 ** -fmt.F))
    err = got - ref
    print(f"  y = {got:.9f}   ln(x) = {ref:.9f}")
    print(f"  abs {err:+.3e}  ({-np.log2(abs(err)):.1f} bits)")


if __name__ == "__main__":
    fmt = DEFAULT
    print(f"LN2 = {fmt.ln2_code}")
    print("coefficients =", [int(c) for c in fmt.c.const(poly.ln1p_coefficients()).code])
    print()

    for ratio in (0.6, 1.05):
        trace(int(round(ratio * 2 ** fmt.F)), fmt)
        print()

    sweep(fmt)
