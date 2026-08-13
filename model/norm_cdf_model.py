"""Golden model for rtl/norm_cdf.sv.

N(x) by Abramowitz-Stegun 26.2.17, with t = 1/(1 + p*|x|) and p = 0.2316419:

    N(|x|) = 1 - phi(|x|) * (b1*t + b2*t^2 + b3*t^3 + b4*t^4 + b5*t^5)
    phi(x) = (1/sqrt(2*pi)) * e^(-x*x/2)

Three rearrangements, all of them free:

1. 1/sqrt(2*pi) is folded into the coefficients at elaboration. That removes a
   multiply, and it drops every coefficient and every internal poly_eval node
   below 1.0 -- max |b'| is 0.7266 against 1.8213 raw -- so the whole polynomial
   needs no integer bit.
2. The symmetry N(-|x|) = 1 - N(|x|) means the negative branch is phi*q' with
   the final subtract skipped. A mux, not a second subtract, and it makes
   N(x) + N(-x) = 1 hold bit-exactly rather than approximately.
3. |x| is clamped to CLAMP before anything else. Past that point N has already
   saturated, so this is exact behaviour and not an approximation, and it is
   also what stops d = 3372 (reachable over the declared input box) from needing
   a format that could hold it.

This composes the recip_unit, poly_eval and exp_unit models rather than
reimplementing them, because the RTL composes those modules. A change to any of
them has to show up here or the mirror has stopped being one.

Specified on absolute error: N is bounded in [0,1] and the pricer adds it to
other terms, so absolute size is what the caller sees -- unlike recip_unit,
which is specified on relative error.

The accuracy ceiling is the formula, not the arithmetic. A&S 26.2.17 is wrong by
up to 7.4517e-08 = 2^-23.68 in exact arithmetic, so spending more than ~23 bits
anywhere in here buys nothing.

Measured error budget, replacing one block at a time with exact arithmetic. The
binding term is not the one you would guess:

    exp_unit at its own defaults (16-bit significand)     17.0 bits
    ... 18 / 20 / 22 fractional bits on the significand   19.2 / 20.3 / 21.1
    ... 24, where it stops paying                         21.4
    recip_unit IDX_W 9 / 10 / 11                          19.3 / 21.4 / 21.6
    t and the polynomial made exact as well                22.5
    the formula alone                                      23.7

exp_unit binds first because its defaults were sized for the discount factor,
where the significand only ever needed 16 bits. Here it sits near 1.0 while N
wants ~24 bits absolute. recip_unit binds second, and IDX_W=10 is free: the seed
ROM was 512 entries in a RAMB18 that holds 1024.

Past ~21.6 bits, t and the polynomial dominate jointly and neither can be fixed
alone -- making one exact is *worse* than leaving both quantised, because their
errors partly cancel (q' is a function of t, so a low t makes a low q'). Buying
the last two bits needs T_F, W_F and U_F all raised together, which pushes three
more signals past the 25-bit DSP port and costs cascades on all of them. The
defaults below stop at 21.1, where the curve flattens.
"""

import math

import exp_unit_model as exp
import poly_eval_model as poly
import recip_unit_model as recip
from fixedpoint import Fmt

P = 0.2316419
B = (0.319381530, -0.356563782, 1.781477937, -1.821255978, 1.330274429)
INV_SQRT_2PI = 1.0 / math.sqrt(2.0 * math.pi)

# b_i / sqrt(2*pi), with the a0 slot poly_eval expects. The series has no
# constant term, which is exactly the a0 = 0 that a fixed six-coefficient
# evaluator already provides -- no need for a degree-4 variant and a spare
# multiply by t.
BFOLD = (0.0,) + tuple(b * INV_SQRT_2PI for b in B)


class NormCdfFormats:
    """The parameters of norm_cdf.sv, and everything derived from them."""

    def __init__(self, IN_W=26, IN_F=22, CLAMP=6.0, P_W=18, P_F=17,
                 T_W=25, T_F=23, C_W=18, C_F=17, W_W=25, W_F=24,
                 U_W=28, U_F=22, PHI_W=32, PHI_F=31, N_W=26, N_F=25,
                 IDX_W=10, EC_W=23, EC_F=21, EW_W=24, EW_F=22):
        self.IN_W, self.IN_F, self.CLAMP = IN_W, IN_F, CLAMP
        self.P_W, self.P_F = P_W, P_F
        self.T_W, self.T_F = T_W, T_F
        self.U_W, self.U_F = U_W, U_F
        self.PHI_W, self.PHI_F = PHI_W, PHI_F
        self.N_W, self.N_F = N_W, N_F

        # |x| keeps the input's binary point, so the clamp is a compare and a
        # sign strip with nothing to round.
        self.AX_W, self.AX_F = IN_W - 1, IN_F
        self.CLAMP_CODE = int(round(CLAMP * 2.0 ** IN_F))
        self.P_CODE = int(round(P * 2.0 ** P_F))

        # t is the reciprocal's output and poly_eval's argument at the same
        # scale, so nothing re-rounds between the two blocks. T_F is one short of
        # T_W-1 because t reaches 1.0 exactly at x = 0 and needs the integer bit.
        self.recip = recip.RecipFormats(Y_W=T_W, Y_F=T_F, IDX_W=IDX_W)
        self.A_F = self.recip.A_F
        self.SHIFT_ARG = P_F + self.AX_F - self.A_F
        self.RND_ARG = 1 << (self.SHIFT_ARG - 1)

        self.c = Fmt(C_W, C_F, "half_up")  # folded coefficients
        self.t_arg = Fmt(T_W, T_F, "half_up")  # t into poly_eval
        self.w = Fmt(W_W, W_F, "half_up")  # poly_eval working format
        self.W_F = W_F

        # u = -|x|^2 / 2. The halving is part of the shift, so it is free.
        self.SHIFT_U = 2 * self.AX_F - U_F + 1
        self.RND_U = 1 << (self.SHIFT_U - 1)

        # exp_unit's MAX_SHIFT is a parameter, not a consequence of U_W/U_F. The
        # format holds -32 because integer bits come in powers of two, but the
        # clamp promises u >= -CLAMP^2/2, so the barrel shifter is built for that
        # and no more. The promise is the clamp's, which is why it is upstream.
        self.MAX_SHIFT = int(math.ceil(0.5 * CLAMP * CLAMP * exp.LOG2E))
        # EC_*/EW_* are exp_unit's own coefficient and working formats. They are
        # separate knobs because exp_unit's defaults were sized for the discount
        # factor, where 16 fractional bits on the significand was ample. Here the
        # significand is near 1 while N wants ~24 bits absolute, so this is the
        # first instance that cares. It is also the term that binds -- see the
        # error budget in the module docstring.
        self.exp = exp.ExpFormats(X_W=U_W, X_F=U_F, L_W=24, L_F=22,
                                  T_W=U_F + 6, T_F=U_F, C_W=EC_W, C_F=EC_F,
                                  W_W=EW_W, W_F=EW_F, Y_W=PHI_W, Y_F=PHI_F)

        self.SHIFT_N = PHI_F + W_F - N_F
        self.RND_N = 1 << (self.SHIFT_N - 1)
        self.ONE_N = 1 << N_F

    @property
    def lo_code(self):
        return -(1 << (self.IN_W - 1))

    @property
    def hi_code(self):
        return (1 << (self.IN_W - 1)) - 1


DEFAULT = NormCdfFormats()


def stages(x, fmt=DEFAULT):
    """Every intermediate for one input code, keyed by the signal in the RTL."""
    x = int(x)
    assert fmt.lo_code <= x <= fmt.hi_code, f"x={x} does not fit in {fmt.IN_W} bits"

    # clock 0. The clamp is the module's only defence: d reaches 3372 over the
    # declared input box and no sane format holds that.
    sign = x < 0
    ax = min(abs(x), fmt.CLAMP_CODE)

    # clock 1-2, branch A. arg = 1 + P*|x|, landing in [1, 2.39] -- two octaves,
    # which is exactly the E_MIN=28, E_MAX=29 window recip_unit defaults to.
    prod_p = fmt.P_CODE * ax
    arg = (1 << fmt.A_F) + ((prod_p + fmt.RND_ARG) >> fmt.SHIFT_ARG)

    # clocks 3-12, the reciprocal. t in [0.418, 1].
    rs = recip.stages(arg, fmt.recip)
    t = rs["y"]

    # clocks 13-19, the polynomial. t is already at the argument format's scale,
    # so const() below recovers this code exactly and adds no rounding.
    q = poly.evaluate_estrin_fixed(BFOLD, t * 2.0 ** -fmt.T_F,
                                   fmt.c, fmt.t_arg, fmt.w)

    # clock 1-2, branch B, five stages later than it could be: the delay is
    # spent on |x| rather than on phi, because |x| is the narrower wire.
    sq = ax * ax
    u = -((sq + fmt.RND_U) >> fmt.SHIFT_U)  # negate last, so the shift is unsigned

    # clocks 8-19, the exponential.
    es = exp.stages(u * 2.0 ** -fmt.U_F, fmt.exp)
    phi = int(es["y"])

    # clock 20-21. N = 1 - phi*q' for x >= 0, and phi*q' for x < 0 -- the
    # symmetry costs a mux, not a subtract. At x = 0 the two agree exactly
    # (phi*q' = 0.5), so the sign of zero never matters.
    prod_n = phi * int(q.code)
    pq = (prod_n + fmt.RND_N) >> fmt.SHIFT_N
    n = pq if sign else fmt.ONE_N - pq

    return {
        "sign": sign, "ax": ax, "prod_p": prod_p, "arg": arg, "t": t,
        "q": int(q.code), "sq": sq, "u": u, "phi": phi,
        "prod_n": prod_n, "pq": pq, "n": n,
    }


def evaluate_fixed(xs, fmt=DEFAULT):
    """The n codes the RTL should produce, for the same x codes."""
    return [stages(x, fmt)["n"] for x in xs]


def evaluate_float(xs, fmt=DEFAULT):
    """Same, as real numbers."""
    return [c * 2.0 ** -fmt.N_F for c in evaluate_fixed(xs, fmt)]


def absolute_error(xs, fmt=DEFAULT):
    """(x, got, error) against the true CDF, clamped where the module clamps."""
    from scipy.stats import norm

    out = []
    for x, got in zip(xs, evaluate_float(xs, fmt)):
        xf = int(x) * 2.0 ** -fmt.IN_F
        ref = float(norm.cdf(max(-fmt.CLAMP, min(fmt.CLAMP, xf))))
        out.append((int(x), got, got - ref))
    return out


def corners(fmt=DEFAULT, span=2):
    """Codes that straddle every boundary the design has.

    Zero, where t reaches 1.0 and the two symmetry branches must agree; the sign
    flip either side of it; the clamp edges, which are the only way to reach the
    saturating path; and 1 + P*|x| = 2, where the reciprocal's exponent steps
    from E_MIN to E_MAX and shift_left changes. That last one is invisible to a
    uniform sweep and is where a wrong window shows up.
    """
    xs = {0}
    lsb = 1
    edges = [fmt.CLAMP_CODE, int(round(2.0 ** fmt.IN_F / P))]  # clamp, arg = 2
    for e in edges:
        for d in range(-span, span + 1):
            xs |= {e + d * lsb, -(e + d * lsb)}
    for d in range(-span, span + 1):
        xs.add(d)
    # a spread over the live range, both signs
    n = 400
    for i in range(n + 1):
        v = int(round(fmt.CLAMP_CODE * i / n))
        xs |= {v, -v}
    # and past the clamp, where every sample must give the same answer
    xs |= {fmt.hi_code, fmt.lo_code, fmt.CLAMP_CODE * 2, -fmt.CLAMP_CODE * 2}
    return sorted(v for v in xs if fmt.lo_code <= v <= fmt.hi_code)


def sweep(fmt=DEFAULT, count=8001):
    """Absolute accuracy over the live range, against scipy."""
    hi = fmt.CLAMP_CODE
    xs = [int(round(-hi + 2 * hi * i / (count - 1))) for i in range(count)]
    res = absolute_error(xs, fmt)
    err = [abs(e) for _, _, e in res]
    worst = max(res, key=lambda r: abs(r[2]))
    bits = lambda v: float("inf") if v <= 0 else -math.log2(v)

    print(f"x in [-{fmt.CLAMP}, {fmt.CLAMP}], {count} codes, clamp {fmt.CLAMP_CODE}")
    print(f"  max abs {max(err):.4e}  ({bits(max(err)):.2f} bits)"
          f"   rms {math.sqrt(sum(v * v for v in err) / len(err)):.4e}")
    print(f"  worst at x = {worst[0] * 2.0 ** -fmt.IN_F:+.9f}"
          f"   got {worst[1]:.12f}")
    print(f"  formula alone is wrong by 7.4517e-08 (23.68 bits), so anything"
          f" past that is the arithmetic")
    print(f"  price/K sees this amplified by S/K <= 30:"
          f" {30 * max(err):.3e} ({bits(30 * max(err)):.2f} bits)")
    return err


def check_symmetry(fmt=DEFAULT, count=2001):
    """N(x) + N(-x) must be exactly 1, by construction rather than by luck.

    x = 0 is excluded because there is only one sample there, not a pair: it
    maps to the positive branch, so the identity would demand N(0) == 1/2
    exactly. It is close but not equal, and that is recip_unit showing through --
    Newton's error rule is eps -> -eps^2, always negative, so 1/1 comes back as
    0.99999976 and never as 1. Reported separately below.
    """
    hi = fmt.CLAMP_CODE
    xs = [int(round(hi * i / (count - 1))) for i in range(count) if i]
    bad = [x for x in xs
           if stages(x, fmt)["n"] + stages(-x, fmt)["n"] != fmt.ONE_N]
    print(f"symmetry N(x) + N(-x) == 1: {len(xs) - len(bad)}/{len(xs)} exact")
    n0 = stages(0, fmt)["n"]
    print(f"  N(0) = {n0 * 2.0 ** -fmt.N_F:.12f}, off 1/2 by "
          f"{n0 - (fmt.ONE_N >> 1)} lsb")
    return bad


def check_monotonic(fmt=DEFAULT, count=20001):
    """A CDF that steps backwards is invisible to an RMS metric and lethal to
    the Phi^-1 round-trip later, so it gets its own check."""
    hi = fmt.CLAMP_CODE
    xs = [int(round(-hi + 2 * hi * i / (count - 1))) for i in range(count)]
    ns = evaluate_fixed(xs, fmt)
    drops = [(xs[i], ns[i], ns[i + 1]) for i in range(len(ns) - 1)
             if ns[i + 1] < ns[i]]
    print(f"monotonic: {len(drops)} backward steps over {len(xs)} codes")
    return drops


def trace(x, fmt=DEFAULT):
    """One code, every intermediate printed. For chasing an RTL mismatch."""
    from scipy.stats import norm

    s = stages(x, fmt)
    print(f"x = {x} = {int(x) * 2.0 ** -fmt.IN_F:+.9f}")
    scale = {"ax": fmt.AX_F, "arg": fmt.A_F, "t": fmt.T_F, "q": fmt.W_F,
             "u": fmt.U_F, "phi": fmt.PHI_F, "pq": fmt.N_F, "n": fmt.N_F}
    for name in ("sign", "ax", "prod_p", "arg", "t", "q", "sq", "u", "phi",
                 "prod_n", "pq", "n"):
        v = int(s[name])
        suffix = f"   = {v * 2.0 ** -scale[name]:+.12f}" if name in scale else ""
        print(f"  {name:9s} {v:>22d}{suffix}")

    got = s["n"] * 2.0 ** -fmt.N_F
    xf = int(x) * 2.0 ** -fmt.IN_F
    ref = float(norm.cdf(max(-fmt.CLAMP, min(fmt.CLAMP, xf))))
    err = got - ref
    print(f"  N = {got:.12f}   true {ref:.12f}")
    print(f"  abs {err:+.3e}"
          f"  ({'exact' if err == 0 else f'{-math.log2(abs(err)):.1f} bits'})")


if __name__ == "__main__":
    fmt = DEFAULT
    print("folded coefficients:")
    for i, b in enumerate(BFOLD):
        print(f"  a{i} = {b:+.12f}   code {fmt.c.const(b).code:>7d}")
    print(f"  max |b'| = {max(abs(b) for b in BFOLD):.9f}, so no integer bit")
    print(f"MAX_SHIFT {fmt.MAX_SHIFT}  SHIFT_ARG {fmt.SHIFT_ARG}  "
          f"SHIFT_U {fmt.SHIFT_U}  SHIFT_N {fmt.SHIFT_N}")
    print()

    trace(0, fmt)
    print()
    trace(int(round(1.0 * 2 ** fmt.IN_F)), fmt)
    print()
    trace(int(round(-2.5 * 2 ** fmt.IN_F)), fmt)
    print()
    sweep(fmt)
    print()
    check_symmetry(fmt)
    check_monotonic(fmt)
