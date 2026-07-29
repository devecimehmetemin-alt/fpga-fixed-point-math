import numpy as np
from scipy.stats import norm

DEGENERATE = 1e-14


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
