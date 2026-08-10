"""Golden model for rtl/sqrt_unit.sv.

sqrt(a) = sqrt(M) * 2^(e/2) for a = M * 2^e out of lzc_norm. e/2 is only an
integer when e is even, so an odd e moves the mantissa instead: M = m for even
e and 2m for odd, which puts M in [1,4) and leaves the shift as floor(e/2).

A midpoint seed table gives 1/sqrt(M) to SEED_BITS, and one Newton step
y1 = y0(1.5 - 0.5*M*y0*y0) squares that error. Factored so the DSP sees it as
C + A*B and the residual stays near zero:

    g0 = M*y0                sqrt(M)
    r  = 0.5*(1 - g0*y0)     near 0, and the 0.5 folds into a shift
    g1 = g0 + g0*r

This mirrors the RTL stage for stage, not a reference sqrt. Specified on
relative error, like recip_unit.
"""

import math


def trunc(v, width):
    """A width cast to an unsigned target: keep the low bits."""
    return v & ((1 << width) - 1)


def wrap(v, width):
    """A width cast to a signed target: keep the low bits, read them signed."""
    v = trunc(v, width)
    return v - (1 << width) if v & (1 << (width - 1)) else v


def seed_table(IDX_W, S_F):
    """Entry i is the inverse square root of slice i's midpoint.

    The top address bit is e[0]: 0 means M in [1,2), 1 means M in [2,4).
    Slicing m rather than M makes the slices geometric in M, so both halves
    carry the same relative error. Must match the RTL's
    $rtoi(1.0/$sqrt(mid) * 2.0**S_F + 0.5), which is round-half-up.
    """
    n, h = 1 << IDX_W, 1 << (IDX_W - 1)
    out = []
    for i in range(n):
        half = i >> (IDX_W - 1)
        mid = (1.0 + ((i & (h - 1)) + 0.5) / h) * (2.0 ** half)
        out.append(math.floor(1.0 / math.sqrt(mid) * 2.0 ** S_F + 0.5))
    return out


class SqrtFormats:
    """The parameters of sqrt_unit.sv, and everything derived from them."""

    def __init__(self, A_W=32, A_F=28, Y_W=24, Y_F=21, IDX_W=10, S_F=17,
                 M_W=25, R_W=18, E_MIN=28, E_MAX=29):
        assert A_F % 2 == 0, "A_F must be even for 2^(-A_F/2) to be a shift"
        self.A_W, self.A_F, self.Y_W, self.Y_F = A_W, A_F, Y_W, Y_F
        self.IDX_W, self.S_F, self.M_W, self.R_W = IDX_W, S_F, M_W, R_W
        self.E_MIN, self.E_MAX = E_MIN, E_MAX

        self.M_F = M_W - 3 # M is in [1,4): two integer bits, one sign
        self.SEED_BITS = IDX_W + 1
        self.NSEED = 1 << IDX_W
        self.HSEED = 1 << (IDX_W - 1)

        self.PROD_W = M_W + S_F + 1
        self.G_F = self.M_F + S_F
        # g_q goes on the wide port: r measures how far g_q*y0 is from 1 and
        # the correction is applied to the full g, so a g_q error passes into
        # the answer undivided. It is not attenuated the way a seed error is.
        self.GB_F = M_W - 3 # g reaches 2.001, so two integer bits
        self.SHIFT_G = self.G_F - self.GB_F

        # as much headroom as the 48-bit accumulator allows on the C port
        self.SHIFT_C = 47 - self.PROD_W
        self.R_F = self.G_F + self.SHIFT_C - self.GB_F
        self.SHIFT_RES = self.GB_F + S_F + 1 - self.R_F # the +1 is the 0.5

        self.G1_F = self.GB_F + self.R_F
        self.G1_W = self.PROD_W + self.SHIFT_C + 1

        self.SHIFT_SPAN = (E_MAX >> 1) - (E_MIN >> 1)
        self.ACC_W = self.G1_W + self.SHIFT_SPAN
        self.SHIFT_RIGHT = self.G1_F - Y_F - ((E_MIN >> 1) - A_F // 2)

        self.PROD2_W = M_W + S_F + 2
        self.ONE = 1 << (self.GB_F + S_F)
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


DEFAULT = SqrtFormats()
SQRT_T = SqrtFormats(Y_W=28, Y_F=25, E_MIN=18, E_MAX=30) # T in [1e-3, 5]


def stages(a, fmt=DEFAULT):
    """Every intermediate for one input code, keyed by the signal in the RTL."""
    a = int(a)
    assert 0 <= a < (1 << fmt.A_W), f"a={a} does not fit in {fmt.A_W} bits"

    # clock 0-1, lzc_norm
    zero = a == 0
    lz = fmt.A_W - a.bit_length() if not zero else fmt.A_W
    m = 0 if zero else a << lz
    e = 0 if zero else fmt.A_W - 1 - lz

    # clock 2. idx takes e[0] and m's fraction: M's leading one moves with the
    # parity, m's never does, so only m gives a uniform slice.
    m_q = m >> (fmt.A_W - (fmt.M_W - 2))
    idx = ((e & 1) << (fmt.IDX_W - 1)) | \
          ((m >> (fmt.A_W - fmt.IDX_W)) & (fmt.HSEED - 1))

    sat_hi = e > fmt.E_MAX
    sat_lo = e < fmt.E_MIN
    if sat_lo:
        shift_left = 0
    elif sat_hi:
        shift_left = fmt.SHIFT_SPAN
    else:
        shift_left = (e >> 1) - (fmt.E_MIN >> 1)

    # clock 3. An odd exponent moves M into [2,4).
    y0 = fmt.SEED[idx]
    M_q = wrap(m_q << (e & 1), fmt.M_W)

    # clock 4-5. DSP1, then g narrowed for the multiplier ports.
    prod = wrap(M_q * y0, fmt.PROD_W)
    g_q = wrap(prod >> fmt.SHIFT_G, fmt.M_W)

    # clock 6-7. DSP2, the residual.
    prod2 = wrap(g_q * y0, fmt.PROD2_W)
    r_raw = wrap(fmt.ONE - prod2, fmt.PROD2_W)

    # clock 8
    r = wrap(r_raw >> fmt.SHIFT_RES, fmt.R_W)

    # clock 9-10. DSP3, g1 = g0 + g0*r, with the full g on the C port.
    prod3 = wrap(g_q * r, fmt.M_W + fmt.R_W)
    g1 = wrap((wrap(prod, fmt.G1_W) << fmt.SHIFT_C) + prod3, fmt.G1_W)

    # clock 11. a = 0 gives e = 0, already below the window, so sqrt(0) = 0
    # falls out of the underflow branch.
    if sat_hi:
        y = fmt.Y_MAX
    elif sat_lo:
        y = 0
    else:
        y = wrap(((wrap(g1, fmt.ACC_W) << shift_left) + fmt.RND_Y)
                 >> fmt.SHIFT_RIGHT, fmt.Y_W)

    return {
        "m": m, "e": e, "zero": zero, "m_q": m_q, "M_q": M_q, "idx": idx,
        "y0": y0, "shift_left": shift_left, "sat_hi": sat_hi, "sat_lo": sat_lo,
        "prod": prod, "g_q": g_q, "prod2": prod2, "r_raw": r_raw, "r": r,
        "prod3": prod3, "g1": g1, "y": y,
    }


def evaluate_fixed(xs, fmt=DEFAULT):
    """The y codes the RTL should produce, for the same a codes."""
    return [stages(x, fmt)["y"] for x in xs]


def evaluate_float(xs, fmt=DEFAULT):
    """Same, as real numbers. Saturated samples come back as None."""
    out = []
    for s in (stages(x, fmt) for x in xs):
        out.append(None if s["sat_hi"] or s["sat_lo"]
                   else s["y"] * 2.0 ** -fmt.Y_F)
    return out


def relative_error(xs, fmt=DEFAULT):
    """(a, got, error) for every sample that is not saturated."""
    out = []
    for x, got in zip(xs, evaluate_float(xs, fmt)):
        if got is None:
            continue
        ref = math.sqrt(int(x) * 2.0 ** -fmt.A_F)
        out.append((int(x), got, (got - ref) / ref))
    return out


def corners(fmt=DEFAULT, per_slice=2):
    """Codes that straddle every boundary the design has.

    The exponent edges, which move shift_left, the parity bit and the
    saturation flags; the seed slice edges, where the table steps and the
    residual reaches its bound; and the codes either side of the declared
    window, which are the only way to reach the saturation paths.
    """
    xs = {0, 1, fmt.lo_code, fmt.hi_code}
    for e in range(fmt.E_MIN - 1, fmt.E_MAX + 2):
        if 0 <= e < fmt.A_W:
            xs |= {1 << e, (1 << e) - 1, (1 << e) + 1}
    for e in range(fmt.E_MIN, fmt.E_MAX + 1):
        # an octave [2^e, 2^(e+1)) split into HSEED slices
        step = 1 << max(0, e - (fmt.IDX_W - 1))
        base = 1 << e
        for i in range(fmt.HSEED):
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
    print(f"  seed alone would give {fmt.SEED_BITS} bits, one Newton step "
          f"reaches {2 * fmt.SEED_BITS - 0.6:.1f}")
    return err


def trace(a, fmt=DEFAULT):
    """One code, every intermediate printed. For chasing an RTL mismatch."""
    s = stages(a, fmt)
    scales = {"m_q": fmt.M_F, "M_q": fmt.M_F, "y0": fmt.S_F, "prod": fmt.G_F,
              "g_q": fmt.GB_F, "r": fmt.R_F, "g1": fmt.G1_F, "y": fmt.Y_F}
    print(f"a = {a} = {a * 2.0 ** -fmt.A_F:.9f}")
    for name in ("m", "e", "zero", "m_q", "M_q", "idx", "y0", "shift_left",
                 "sat_hi", "sat_lo", "prod", "g_q", "prod2", "r_raw", "r",
                 "prod3", "g1", "y"):
        v = int(s[name])
        suffix = f"   = {v * 2.0 ** -scales[name]:+.9f}" if name in scales else ""
        print(f"  {name:11s} {v:>18d}{suffix}")

    if not s["sat_hi"] and not s["sat_lo"]:
        got = s["y"] * 2.0 ** -fmt.Y_F
        ref = math.sqrt(a * 2.0 ** -fmt.A_F)
        rel = (got - ref) / ref
        print(f"  sqrt(a) = {got:.9f}   true {ref:.9f}")
        print(f"  rel {rel:+.3e}  ({-math.log2(abs(rel)):.1f} bits)")


if __name__ == "__main__":
    fmt = DEFAULT
    print(f"SEED[0] {fmt.SEED[0]}  SEED[{fmt.HSEED}] {fmt.SEED[fmt.HSEED]}  "
          f"SEED[{fmt.NSEED-1}] {fmt.SEED[-1]}  max {max(fmt.SEED)} "
          f"< 2^{fmt.S_F} = {1 << fmt.S_F}")
    print(f"M_F {fmt.M_F}  G_F {fmt.G_F}  GB_F {fmt.GB_F}  SHIFT_G {fmt.SHIFT_G}  "
          f"R_F {fmt.R_F}  SHIFT_RES {fmt.SHIFT_RES}")
    print(f"G1_F {fmt.G1_F}  SHIFT_C {fmt.SHIFT_C}  SHIFT_RIGHT {fmt.SHIFT_RIGHT}  "
          f"SPAN {fmt.SHIFT_SPAN}")
    print()

    trace(round(2.0 * 2 ** fmt.A_F), fmt)
    print()
    sweep(fmt)
    print()
    sweep(SQRT_T)
