"""Golden model for rtl/bs_pricer_top.sv -- the fixed-point pricer.

    v   = sigma*sqrt(T)
    d1  = [ln(S/K) + (r + sigma^2/2)*T] / v
    d2  = d1 - v
    C   = S*N(d1) - K*e^(-rT)*N(d2)
    P   = K*e^(-rT)*N(-d2) - S*N(-d1)

bs_pricer_model.py is the float reference; this is the spec the RTL must match
bit for bit. It composes the unit models rather than reimplementing them,
because the RTL composes those modules.

Four structural decisions, all visible below:

1. K-normalised interior. S and K enter, but the datapath carries x = S/K and
   produces C/K, multiplying by K once at the very end. x has to be computed
   anyway for ln(S/K), so this costs nothing and it deletes S's delay line --
   S would otherwise have to survive ~60 clocks to reach the final combine.

2. d1 and d2 saturate to CLAMP *in the pricer*, not in norm_cdf. Measured over
   the declared box |d| reaches 10756, and norm_cdf's input format holds 8, so
   the value would wrap on the way in and the clamp inside norm_cdf would never
   see the number it was meant to catch. norm_cdf's own clamp becomes a no-op
   here, which is the correct relationship: the module guarantees it for any
   caller, this caller has already guaranteed it.

3. d2 is formed from the *unsaturated* d1. Saturating first and subtracting
   after is a different function wherever d1 is past the clamp, and it is the
   version you write by accident. bs_pricer_model.d1_d2 clips both after the
   subtraction for exactly this reason.

4. No degenerate bypass. The spec calls for one as v -> 0, but over RANGES
   v >= 3.16e-04, and recip_unit saturating below its exponent window drives
   |d| up rather than to a nan -- so the clamp already delivers discounted
   intrinsic. The bypass is needed only if the declared box is widened to
   include sigma = 0 or T = 0.

Measured: 19.3 bits on C/K over the box. The binding term is x = S/K, and the
isolation is unambiguous -- making lnx, 1/v, hvT and v all exact simultaneously
moves the error by nothing at all, while x alone accounts for it. inv_k tops out
at 23.0 bits relative and stays there whether IK_F is 32 or 35, whether the seed
is 2^12 or 2^14, and whether X_F is 30 or 34, so the limit is upstream of
recip_unit: K's own quantisation at SK_F. x then carries that 23 bits into
C/K = x*N(d1) - ..., where it is amplified by x itself, which reaches 29 over
RANGES. 23 - log2(29) = 18.1, against 19.3 measured.

So the lever is SK_F, the input format, and nothing inside the datapath. That is
worth knowing before anyone widens a working format hoping to buy accuracy.

Two hypotheses that measurement killed, recorded so they are not retried:
conditioning (d1 = num/v amplifies num's error by 1/v <= 3162) does not bind,
because wherever 1/v is large d1 is past the clamp and N has saturated; and the
A&S formula ceiling does not bind either -- swapping the reference from erf to
A&S-in-float64 changes the error by 10%, not by orders of magnitude.
"""

import math

import exp_unit_model as exp
import ln_unit_model as ln
import norm_cdf_model as nc
import recip_unit_model as recip
import sqrt_unit_model as sqrt


def _sat(code, width):
    """Saturating cast to a signed width. Wrapping here would be silent."""
    hi = (1 << (width - 1)) - 1
    lo = -(1 << (width - 1))
    return hi if code > hi else (lo if code < lo else code)


def _rshift(code, shift):
    """Round-half-up right shift, the same one every unit uses."""
    if shift <= 0:
        return code << -shift
    return (code + (1 << (shift - 1))) >> shift


class PricerFormats:
    """The parameters of bs_pricer_top.sv, and everything derived from them."""

    def __init__(self, SK_W=32, SK_F=22, T_W=32, T_F=28, SG_W=24, SG_F=22,
                 R_W=24, R_F=22, IK_W=32, IK_F=32, X_W=36, X_F=30,
                 LX_W=28, LX_F=24, ST_W=28, ST_F=25, V_W=28, V_F=25,
                 IV_W=32, IV_F=19, HV_W=26, HV_F=24, HT_W=28, HT_F=24,
                 NM_W=29, NM_F=24, D_W=40, IN_W=26, IN_F=22,
                 U_W=28, U_F=22, DC_W=32, DC_F=31, CN_W=32, CN_F=26,
                 PN_W=28, PN_F=26, OUT_W=40, OUT_F=22,
                 K_IDX=12, V_IDX=10, EC_W=23, EC_F=21, EW_W=24, EW_F=22):
        for k, v in list(locals().items()):
            if k != "self":
                setattr(self, k, v)

        # 1/K. K in [10,300] at SK_F=22 puts the exponent in [25,30], six
        # octaves. IK_F=32 rather than the width minus the integer bits: 1/K
        # never reaches 1/8, so the format spends every bit on the fraction and
        # still holds the top of the range.
        self.recip_k = recip.RecipFormats(A_W=SK_W, A_F=SK_F, Y_W=IK_W,
                                          Y_F=IK_F, IDX_W=K_IDX,
                                          E_MIN=25, E_MAX=30)
        # K_IDX=12 rather than the 10 used elsewhere. This seed is the only one
        # whose relative error lands on ln, and d(ln x) = dx/x carries it
        # through undiminished, straight into the 1/v amplifier below.
        self.SHIFT_X = SK_F + IK_F - X_F

        # ln(S/K). X_F=30 on a value spanning ten octaves is the compromise:
        # at x = 1/30 an lsb is 2^-25.1 relative, which is what caps this path.
        self.ln = ln.LnFormats(W=X_W, F=X_F, X_W=25, X_F=24, C_W=25, C_F=24,
                               W_W=25, W_F=LX_F, L_W=26, Y_W=LX_W)

        # sqrt(T). T in [1e-3,5] at T_F=28 puts the exponent in [18,30], and
        # A_F must be even for 2^(-A_F/2) to stay a shift.
        self.sqrt_t = sqrt.SqrtFormats(A_W=T_W, A_F=T_F, Y_W=ST_W, Y_F=ST_F,
                                       IDX_W=10, E_MIN=18, E_MAX=30)
        self.SHIFT_V = SG_F + ST_F - V_F

        # 1/v. v in [3.16e-4, 3.354] at V_F=25 is fourteen octaves -- the widest
        # window any unit in the project is asked for, and the reason recip's
        # output shifter is the deepest instance here.
        self.recip_v = recip.RecipFormats(A_W=IV_W, A_F=V_F, Y_W=IV_W,
                                          Y_F=IV_F, IDX_W=V_IDX,
                                          E_MIN=13, E_MAX=26)

        self.SHIFT_HV = 2 * SG_F + 1 - HV_F  # sigma^2/2: the halving is free
        self.SHIFT_HT = HV_F + T_F - HT_F
        self.SHIFT_NM = LX_F - NM_F  # lnx and hvT meet at NM_F
        self.SHIFT_D = NM_F + IV_F - IN_F
        self.SHIFT_V2D = V_F - IN_F  # v down to d's scale for d2 = d1 - v

        # -rT, and e^(-rT). u in [-0.5, 0.1] means the exponent of the result is
        # 0 or -1, so MAX_SHIFT is 1 -- exp_unit's default, and the one place in
        # the project where its defaults were already right.
        self.SHIFT_U = R_F + T_F - U_F
        self.exp = exp.ExpFormats(X_W=U_W, X_F=U_F, L_W=24, L_F=22,
                                  T_W=U_F + 6, T_F=U_F, C_W=EC_W, C_F=EC_F,
                                  W_W=EW_W, W_F=EW_F, Y_W=DC_W, Y_F=DC_F)
        self.MAX_SHIFT = 1

        self.ncdf = nc.NormCdfFormats(IN_W=IN_W, IN_F=IN_F, EC_W=EC_W,
                                      EC_F=EC_F, EW_W=EW_W, EW_F=EW_F)
        self.N_F = self.ncdf.N_F
        self.ONE_N = self.ncdf.ONE_N
        self.CLAMP_CODE = self.ncdf.CLAMP_CODE

        self.SHIFT_XN = X_F + self.N_F - CN_F  # x*N(d1)
        self.SHIFT_DN = DC_F + self.N_F - CN_F  # disc*N(d2)
        self.SHIFT_PN = DC_F + self.N_F - PN_F
        self.SHIFT_XP = X_F + self.N_F - PN_F
        self.SHIFT_C = SK_F + CN_F - OUT_F  # K * (C/K)
        self.SHIFT_P = SK_F + PN_F - OUT_F

    def code(self, value, frac):
        return int(math.floor(value * (2.0 ** frac) + 0.5))


DEFAULT = PricerFormats()


def stages(S, K, r, sigma, T, fmt=DEFAULT):
    """Every intermediate for one input, keyed by the signal in the RTL.

    S, K, r, sigma, T are integer codes at SK_F, SK_F, R_F, SG_F, T_F.
    """
    S, K, r, sigma, T = int(S), int(K), int(r), int(sigma), int(T)

    # clocks 0-10, branch K: 1/K, then x = S/K.
    inv_k = recip.stages(K, fmt.recip_k)["y"]
    x = _rshift(S * inv_k, fmt.SHIFT_X)

    # clocks 11-21, branch K continued: ln(x). ln_unit reads x as Q(W-F).F and
    # subtracts F from the exponent itself, so x goes in as a raw code.
    lnx = int(ln.stages([x], fmt.ln)["y"][0])

    # clocks 0-11, branch T: sqrt(T), then v = sigma*sqrt(T).
    sq_t = sqrt.stages(T, fmt.sqrt_t)["y"]
    v = _rshift(sigma * sq_t, fmt.SHIFT_V)

    # clocks 12-21, branch T continued: 1/v.
    inv_v = recip.stages(v, fmt.recip_v)["y"]

    # clocks 0-2, branch R: (r + sigma^2/2)*T. Cheap and early, so it waits.
    hv = (r << (fmt.HV_F - fmt.R_F)) + _rshift(sigma * sigma, fmt.SHIFT_HV)
    hv_t = _rshift(hv * T, fmt.SHIFT_HT)

    # clock 22: the numerator. lnx and hv_t are both at LX_F/HT_F = NM_F.
    num = _rshift(lnx, fmt.SHIFT_NM) + hv_t

    # clock 23-24: d1 wide, then d2 from the *unsaturated* d1, then both
    # saturated to the CDF's input format at the clamp.
    d1_wide = _rshift(num * inv_v, fmt.SHIFT_D)
    d2_wide = d1_wide - _rshift(v, fmt.SHIFT_V2D)
    clip = lambda c: max(-fmt.CLAMP_CODE, min(fmt.CLAMP_CODE, c))
    d1 = clip(_sat(d1_wide, fmt.D_W))
    d2 = clip(_sat(d2_wide, fmt.D_W))

    # clocks 25-46: the two CDFs, and in parallel the discount factor.
    nd1 = nc.stages(d1, fmt.ncdf)["n"]
    nd2 = nc.stages(d2, fmt.ncdf)["n"]

    u = -_rshift(r * T, fmt.SHIFT_U)
    disc = int(exp.stages(u * 2.0 ** -fmt.U_F, fmt.exp)["y"])

    # clocks 47-49: combine, K-normalised, then scale by K once.
    # 1 - N is exact: the CDF's output format has 1.0 in it by construction.
    call_n = (_rshift(x * nd1, fmt.SHIFT_XN)
              - _rshift(disc * nd2, fmt.SHIFT_DN))
    put_n = (_rshift(disc * (fmt.ONE_N - nd2), fmt.SHIFT_PN)
             - _rshift(x * (fmt.ONE_N - nd1), fmt.SHIFT_XP))

    call_raw = _rshift(K * call_n, fmt.SHIFT_C)
    put_raw = _rshift(K * put_n, fmt.SHIFT_P)

    # Deep out of the money both terms are nearly equal and the difference lands
    # a few lsb below zero. A negative option price is not a rounding detail to
    # the caller, so the output floors at zero -- a sign test and a mux. Parity
    # is asserted on call_raw/put_raw, before this, because the floor is exactly
    # where C - P would stop equalling S - K*disc.
    call = max(call_raw, 0)
    put = max(put_raw, 0)

    return {
        "inv_k": inv_k, "x": x, "lnx": lnx, "sq_t": sq_t, "v": v,
        "inv_v": inv_v, "hv": hv, "hv_t": hv_t, "num": num,
        "d1_wide": d1_wide, "d2_wide": d2_wide, "d1": d1, "d2": d2,
        "nd1": nd1, "nd2": nd2, "u": u, "disc": disc,
        "call_n": call_n, "put_n": put_n,
        "call_raw": call_raw, "put_raw": put_raw, "call": call, "put": put,
    }


def encode(S, K, r, sigma, T, fmt=DEFAULT):
    """Reals to the input codes, in the order the RTL takes them."""
    return (fmt.code(S, fmt.SK_F), fmt.code(K, fmt.SK_F), fmt.code(r, fmt.R_F),
            fmt.code(sigma, fmt.SG_F), fmt.code(T, fmt.T_F))


def evaluate_fixed(inputs, fmt=DEFAULT):
    """The (call, put) codes the RTL should produce, for the same input codes."""
    out = []
    for row in inputs:
        s = stages(*row, fmt=fmt)
        out.append((s["call"], s["put"]))
    return out


def evaluate_float(inputs, fmt=DEFAULT):
    return [(c * 2.0 ** -fmt.OUT_F, p * 2.0 ** -fmt.OUT_F)
            for c, p in evaluate_fixed(inputs, fmt)]


def parity_residual(inputs, fmt=DEFAULT):
    """C - P - (S - K*e^(-rT)), on the module's own numbers.

    Independent of any reference: it is a statement the pricer makes about
    itself, so it holds on the RTL even where the price is inaccurate.
    """
    out = []
    for row in inputs:
        s = stages(*row, fmt=fmt)
        S = row[0] * 2.0 ** -fmt.SK_F
        K = row[1] * 2.0 ** -fmt.SK_F
        disc = s["disc"] * 2.0 ** -fmt.DC_F
        c = s["call_raw"] * 2.0 ** -fmt.OUT_F
        p = s["put_raw"] * 2.0 ** -fmt.OUT_F
        out.append(c - p - (S - K * disc))
    return out
