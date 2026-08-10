import numpy as np
from scipy.stats import norm

DEGENERATE = 1e-14

# The declared operating envelope. Every study and every testbench draws from
# here and nowhere else: the clamp bound below, the exponent range norm_cdf
# hands exp_unit, and the fixed-point formats downstream are all derived from
# these five intervals, so a literal copied elsewhere is a silent spec fork.
#
# The corner sigma=0.01 with T=1e-3 is deliberate, not an oversight. It drives
# sigma*sqrt(T) to 3e-4 and |d| past 3000, which is what exercises the
# saturating path through recip_unit. Removing it would hide the case the
# hardware most needs to survive.
RANGES = {
    "S": (10.0, 300.0),  # spot
    "K": (10.0, 300.0),  # strike
    "r": (-0.02, 0.10),  # risk-free rate, negative rates included
    "sigma": (0.01, 1.50),  # volatility, index through distressed
    "T": (1e-3, 5.0),  # maturity in years, ~8 hours to 5 years
}

# |d1|, |d2| are clamped to this before the CDF. Not a range check -- N() has
# already saturated by here, so the clamp is exact behaviour rather than an
# approximation, and it delivers the sigma*sqrt(T) -> 0 degenerate case for
# free: both CDFs go to 1 and the price falls out as S - K*exp(-r*T).
#
# 6.0 is chosen by measurement, not by N(-6). The error the clamp injects is
# (S/K)*(1 - N(D_CLAMP)), and S/K reaches 30 over RANGES, so the price error is
# ~30x the tail probability. Measured over the box: 2.8e-8 relative to K, an
# eighth of an lsb at Q.22. At 5.657 -- which would land the exponent argument
# exactly on exp_unit's characterised [-16, 0] -- it is 2.2e-7, a full lsb.
D_CLAMP = 6.0


ARGS = ("S", "K", "r", "sigma", "T")


def sample(n, rng):
    """n draws from RANGES, in the argument order the pricer takes."""
    return tuple(rng.uniform(*RANGES[k], n) for k in ARGS)


def d1_d2(S, K, r, sigma, T):
    v = sigma * np.sqrt(T)
    with np.errstate(divide="ignore", invalid="ignore"):
        d1 = (np.log(S / K) + (r + 0.5 * sigma * sigma) * T) / v
    return d1, d1 - v, v


def call(S, K, r, sigma, T):
    d1, d2, v = d1_d2(S, K, r, sigma, T)
    disc_K = K * np.exp(-r * T)
    with np.errstate(invalid="ignore"):
        px = S * norm.cdf(d1) - disc_K * norm.cdf(d2)
    return np.where(v > DEGENERATE, px, np.maximum(S - disc_K, 0.0))[()]


def put(S, K, r, sigma, T):
    d1, d2, v = d1_d2(S, K, r, sigma, T)
    disc_K = K * np.exp(-r * T)
    with np.errstate(invalid="ignore"):
        px = disc_K * norm.cdf(-d2) - S * norm.cdf(-d1)
    return np.where(v > DEGENERATE, px, np.maximum(disc_K - S, 0.0))[()]


def greeks(S, K, r, sigma, T, kind="call"):
    d1, d2, v = d1_d2(S, K, r, sigma, T)
    disc_K = K * np.exp(-r * T)
    live = v > DEGENERATE

    with np.errstate(invalid="ignore", divide="ignore"):
        pdf_d1 = norm.pdf(d1)
        decay = np.where(live, -S * pdf_d1 * sigma / (2.0 * np.sqrt(T)), 0.0)
        out = {
            "gamma": np.where(live, pdf_d1 / (S * v), 0.0),
            "vega": np.where(live, S * pdf_d1 * np.sqrt(T), 0.0),
        }
        if kind == "call":
            itm = np.where(S > disc_K, 1.0, 0.0)
            out["delta"] = np.where(live, norm.cdf(d1), itm)
            out["theta"] = decay - r * disc_K * np.where(live, norm.cdf(d2), itm)
            out["rho"] = np.where(live, K * T * np.exp(-r * T) * norm.cdf(d2), 0.0)
        else:
            itm = np.where(S < disc_K, 1.0, 0.0)
            out["delta"] = np.where(live, norm.cdf(d1) - 1.0, -itm)
            out["theta"] = decay + r * disc_K * np.where(live, norm.cdf(-d2), itm)
            out["rho"] = np.where(live, -K * T * np.exp(-r * T) * norm.cdf(-d2), 0.0)

    return {k: x[()] for k, x in out.items()}


def parity_residual(S, K, r, sigma, T):
    return call(S, K, r, sigma, T) - put(S, K, r, sigma, T) - (S - K * np.exp(-r * T))


def call_normalised(S_over_K, r, sigma, T):
    return call(S_over_K, 1.0, r, sigma, T)
