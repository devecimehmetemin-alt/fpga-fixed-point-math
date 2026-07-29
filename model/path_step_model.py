import numpy as np

from fixedpoint import VARS


def step_float(S, dW, rh, sh2):
    mult1 = sh2 * dW
    sum1 = rh + mult1
    mult2 = S * sum1
    return S + mult2


def step_fixed(S, dW, rh, sh2, F, coarse=False):
    mult1 = F.mult1.cast(sh2 * dW)
    sum1 = F.sum1.cast(rh + mult1)
    if coarse:
        return F.Sc.cast(S + F.mult2c.cast(S * sum1))
    return F.S.cast(S + F.mult2.cast(S * sum1))


def measure_ranges(S0, r, sigma, T, M, n_paths, rng, pad=1.25):
    h = T / M
    rh, sh2 = r * h, sigma * np.sqrt(h)
    peak = {v: 0.0 for v in VARS}
    peak["rh"], peak["sh2"] = abs(rh), abs(sh2)

    S = np.full(n_paths, S0)
    Sc = np.full(n_paths, S0)
    for _ in range(M // 2):
        for z in (rng.standard_normal(n_paths), rng.standard_normal(n_paths)):
            m1, s1 = sh2 * z, rh + sh2 * z
            m2 = S * s1
            S = S + m2
            for k, arr in (("dW", z), ("mult1", m1), ("sum1", s1), ("mult2", m2), ("S", S)):
                peak[k] = max(peak[k], float(np.abs(arr).max()))
        zc = rng.standard_normal(n_paths) * np.sqrt(2.0)
        m2c = Sc * (2.0 * rh + sh2 * zc)
        Sc = Sc + m2c
        peak["mult2c"] = max(peak["mult2c"], float(np.abs(m2c).max()))
        peak["Sc"] = max(peak["Sc"], float(np.abs(Sc).max()))

    return {k: v * pad for k, v in peak.items()}
