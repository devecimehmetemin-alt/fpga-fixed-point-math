"""Golden model for rtl/recip_unit.sv.

1/a = (1/m) * 2^-e for a = m * 2^e out of lzc_norm. A midpoint seed table gives
1/m to SEED_BITS, and one Newton step y1 = y0(2 - m*y0) squares that error, so
the answer lands at 2*SEED_BITS bits.

The step is factored as y1 = y0 + y0*r with r = 1 - m*y0, since
y0(2 - m*y0) = y0(1 + (1 - m*y0)) = y0(1 + r). Same value, different bits: r is
near 0 and left-aligns, while 2 - m*y0 is near 1 and hides the correction in its
low bits.

This mirrors the RTL stage for stage, not a reference reciprocal. Everything is
plain Python ints because that is what the hardware holds -- the binary points
live in the format constants, never in a value.

Specified on relative error, unlike ln_unit: 1/a spans orders of magnitude and
the caller cares about significant digits, not absolute size.
"""

import math


def trunc(v, width):
    """A width cast to an unsigned target: keep the low bits."""
    return v & ((1 << width) - 1)


def wrap(v, width):
    """A width cast to a signed target: keep the low bits, read them signed.

    The RTL never intends to wrap here. Modelling it anyway means a mismatch is
    the design's, not a Python int that quietly grew wider than the wire.
    """
    v = trunc(v, width)
    return v - (1 << width) if v & (1 << (width - 1)) else v


def seed_table(IDX_W, S_F):
    """Entry i is the reciprocal of slice i's midpoint.

    Must match the RTL's $rtoi(1.0/mid * 2.0**S_F + 0.5) exactly. $rtoi
    truncates toward zero and the value is positive, so that is round-half-up --
    not Python's round(), which is half-to-even.
    """
    n = 1 << IDX_W
    return [math.floor(1.0 / (1.0 + (i + 0.5) / n) * 2.0 ** S_F + 0.5)
            for i in range(n)]


class RecipFormats:
    """The parameters of recip_unit.sv, and everything derived from them."""

    def __init__(self, A_W=32, A_F=28, Y_W=24, Y_F=22, IDX_W=9, S_F=17,
                 M_W=25, R_W=18, E_MIN=28, E_MAX=29):
        self.A_W, self.A_F, self.Y_W, self.Y_F = A_W, A_F, Y_W, Y_F
        self.IDX_W, self.S_F, self.M_W, self.R_W = IDX_W, S_F, M_W, R_W
        self.E_MIN, self.E_MAX = E_MIN, E_MAX

        self.M_F = M_W - 2 # one integer bit, one sign
        self.SEED_BITS = IDX_W + 1
        self.R_F = self.SEED_BITS + R_W - 2
        self.SHIFT_RES = self.M_F + S_F - self.R_F
        self.NSEED = 1 << IDX_W
        self.PROD_W = M_W + S_F + 1
        self.PROD2_W = S_F + R_W + 1
        self.Y1_W = S_F + self.R_F + 2
        self.ACC_W = self.Y1_W + E_MAX - E_MIN
        self.SHIFT_FIX = S_F + self.R_F - A_F - Y_F + E_MIN
        self.SHIFT_RIGHT = self.SHIFT_FIX + E_MAX - E_MIN

        self.ONE = 1 << (self.M_F + S_F)
        self.RND_Y = 1 << (self.SHIFT_RIGHT - 1)
        self.Y_MAX = (1 << (Y_W - 1)) - 1
        self.SEED = seed_table(IDX_W, S_F)

    @property
    def lo_code(self):
        """Smallest a inside the declared exponent window."""
        return 1 << self.E_MIN

    @property
    def hi_code(self):
        """Largest a inside it."""
        return (1 << (self.E_MAX + 1)) - 1


DEFAULT = RecipFormats() # the norm_cdf instance, a = 1 + p*x
SIGMA_SQRT_T = RecipFormats(Y_W=32, Y_F=21, E_MIN=18, E_MAX=29)


def stages(a, fmt=DEFAULT):
    """Every intermediate for one input code, keyed by the signal in the RTL."""
    a = int(a)
    assert 0 <= a < (1 << fmt.A_W), f"a={a} does not fit in {fmt.A_W} bits"

    # clock 0-1, lzc_norm
    zero = a == 0
    lz = fmt.A_W - a.bit_length() if not zero else fmt.A_W
    m = 0 if zero else a << lz
    e = 0 if zero else fmt.A_W - 1 - lz

    # clock 2. m_q keeps the leading one because the multiplier needs the value;
    # idx skips it because a constant bit is a dead address line.
    m_q = m >> (fmt.A_W - (fmt.M_W - 1))
    idx = (m >> (fmt.A_W - 1 - fmt.IDX_W)) & (fmt.NSEED - 1)

    sat_hi = e < fmt.E_MIN
    sat_lo = e > fmt.E_MAX
    if sat_hi:
        shift_left = fmt.E_MAX - fmt.E_MIN
    elif sat_lo:
        shift_left = 0
    else:
        shift_left = fmt.E_MAX - e

    # clock 3
    y0 = fmt.SEED[idx]

    # clock 4-5. One DSP: the product in MREG, ONE - prod in the post-adder.
    prod = wrap(m_q * y0, fmt.PROD_W)
    r_raw = wrap(fmt.ONE - prod, fmt.PROD_W)

    # clock 6. r_raw's top SEED_BITS bits are sign, so the slice starts below
    # them. An arithmetic shift right is floor, which is what Python >> does.
    r = wrap(r_raw >> fmt.SHIFT_RES, fmt.R_W)

    # clock 7-8. The second DSP, y1 = y0 + y0*r.
    prod2 = wrap(y0 * r, fmt.PROD2_W)
    y1 = wrap((wrap(y0, fmt.Y1_W) << fmt.R_F) + prod2, fmt.Y1_W)

    # clock 9. Variable left shift then one fixed right shift, so the rounding
    # constant sits at a single known position.
    if zero or sat_hi:
        y, ovf = fmt.Y_MAX, 1
    elif sat_lo:
        y, ovf = 0, 0
    else:
        y = wrap(((wrap(y1, fmt.ACC_W) << shift_left) + fmt.RND_Y) >> fmt.SHIFT_RIGHT,
                 fmt.Y_W)
        ovf = 0

    return {
        "m": m, "e": e, "zero": zero, "m_q": m_q, "idx": idx, "y0": y0,
        "shift_left": shift_left, "sat_hi": sat_hi, "sat_lo": sat_lo,
        "prod": prod, "r_raw": r_raw, "r": r, "prod2": prod2, "y1": y1,
        "y": y, "ovf": ovf,
    }


def evaluate_fixed(xs, fmt=DEFAULT):
    """The (y, ovf) pairs the RTL should produce, for the same a codes."""
    return [(s["y"], s["ovf"]) for s in (stages(x, fmt) for x in xs)]


def evaluate_float(xs, fmt=DEFAULT):
    """Same, as real numbers. Saturated samples come back as None."""
    out = []
    for s in (stages(x, fmt) for x in xs):
        out.append(None if s["ovf"] or s["sat_lo"] else s["y"] * 2.0 ** -fmt.Y_F)
    return out


def relative_error(xs, fmt=DEFAULT):
    """(a, got, error) for every sample that is not saturated."""
    out = []
    for x, got in zip(xs, evaluate_float(xs, fmt)):
        if got is None:
            continue
        ref = 2.0 ** fmt.A_F / int(x)
        out.append((int(x), got, (got - ref) / ref))
    return out


def corners(fmt=DEFAULT, per_slice=2):
    """Codes that straddle every boundary the design has.

    The exponent edges, which move shift_left and the saturation flags; the seed
    slice edges, where the table steps and the residual reaches its bound; and
    the codes either side of the declared window, which are the only way to
    reach the saturation paths at all.
    """
    xs = {0, 1, fmt.lo_code, fmt.hi_code}
    for e in range(fmt.E_MIN - 1, fmt.E_MAX + 2):
        if 0 <= e < fmt.A_W:
            xs |= {1 << e, (1 << e) - 1, (1 << e) + 1}
    for e in range(fmt.E_MIN, fmt.E_MAX + 1):
        # an octave [2^e, 2^(e+1)) split into NSEED slices
        step = 1 << max(0, e - fmt.IDX_W)
        base = 1 << e
        for i in range(fmt.NSEED):
            for d in range(per_slice):
                v = base + i * step + d
                if v < (1 << fmt.A_W):
                    xs.add(v)
    return sorted(v for v in xs if 0 <= v < (1 << fmt.A_W))


def sweep(fmt=DEFAULT, count=20001):
    """Relative accuracy across the declared window."""
    lo, hi = fmt.lo_code, fmt.hi_code
    step = max(1, (hi - lo) // count)
    xs = list(range(lo, hi + 1, step))
    res = relative_error(xs, fmt)
    err = [abs(r) for _, _, r in res]
    worst = max(res, key=lambda t: abs(t[2]))

    print(f"a in [{lo * 2.0 ** -fmt.A_F:.6g}, {hi * 2.0 ** -fmt.A_F:.6g}], "
          f"{len(xs)} codes, e in [{fmt.E_MIN}, {fmt.E_MAX}]")
    print(f"  max rel {max(err):.3e}  ({-math.log2(max(err)):.1f} bits)"
          f"   rms {math.sqrt(sum(v * v for v in err) / len(err)):.3e}")
    print(f"  worst at a = {worst[0]} = {worst[0] * 2.0 ** -fmt.A_F:.9f}")
    print(f"  seed alone would give {fmt.SEED_BITS} bits, "
          f"one Newton step doubles it to {2 * fmt.SEED_BITS}")
    return err


def trace(a, fmt=DEFAULT):
    """One code, every intermediate printed. For chasing an RTL mismatch."""
    s = stages(a, fmt)
    print(f"a = {a} = {a * 2.0 ** -fmt.A_F:.9f}")
    for name in ("m", "e", "zero", "m_q", "idx", "y0", "shift_left", "sat_hi",
                 "sat_lo", "prod", "r_raw", "r", "prod2", "y1", "y", "ovf"):
        v = int(s[name])
        suffix = ""
        if name == "m_q":
            suffix = f"   = {v * 2.0 ** -fmt.M_F:+.9f}"
        elif name == "y0":
            suffix = f"   = {v * 2.0 ** -fmt.S_F:+.9f}"
        elif name == "r":
            suffix = f"   = {v * 2.0 ** -fmt.R_F:+.3e}"
        elif name == "y1":
            suffix = f"   = {v * 2.0 ** -(fmt.S_F + fmt.R_F):+.9f}"
        elif name == "y":
            suffix = f"   = {v * 2.0 ** -fmt.Y_F:+.9f}"
        print(f"  {name:10s} {v:>16d}{suffix}")

    if not s["ovf"] and not s["sat_lo"]:
        got = s["y"] * 2.0 ** -fmt.Y_F
        ref = 2.0 ** fmt.A_F / a
        rel = (got - ref) / ref
        print(f"  1/a = {got:.9f}   true {ref:.9f}")
        print(f"  rel {rel:+.3e}  ({-math.log2(abs(rel)):.1f} bits)")


if __name__ == "__main__":
    fmt = DEFAULT
    print(f"SEED[0] {fmt.SEED[0]}  SEED[{fmt.NSEED-1}] {fmt.SEED[-1]}  "
          f"max {max(fmt.SEED)} < 2^{fmt.S_F} = {1 << fmt.S_F}")
    print(f"R_F {fmt.R_F}  SHIFT_RES {fmt.SHIFT_RES}  SHIFT_FIX {fmt.SHIFT_FIX}  "
          f"SHIFT_RIGHT {fmt.SHIFT_RIGHT}  ONE 1<<{fmt.M_F + fmt.S_F}")
    print()

    trace(round(1.2316419 * 2 ** fmt.A_F), fmt)
    print()
    sweep(fmt)
    print()
    sweep(SIGMA_SQRT_T)
