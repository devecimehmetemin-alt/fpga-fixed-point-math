import numpy as np

from path_step_model import step_fixed, step_float


def simulate(S0, r, sigma, T, M, n_paths, rng):
    h = T / M
    rh = r * h
    # dW is a unit normal, not the Brownian increment
    sh2 = sigma * np.sqrt(h)

    S = np.full(n_paths, S0, dtype=float)
    for _ in range(M):
        S = step_float(S, rng.standard_normal(n_paths), rh, sh2)
    return S


def simulate_coupled(S0, r, sigma, T, M_fine, n_paths, rng):
    if M_fine % 2:
        raise ValueError("M_fine must be even to pair increments")

    h = T / M_fine
    rh = r * h
    sh2 = sigma * np.sqrt(h)

    S = np.full(n_paths, S0, dtype=float)
    Sc = np.full(n_paths, S0, dtype=float)

    for _ in range(M_fine // 2):
        z1 = rng.standard_normal(n_paths)
        z2 = rng.standard_normal(n_paths)

        S = step_float(S, z1, rh, sh2)
        S = step_float(S, z2, rh, sh2)

        Sc = step_float(Sc, z1 + z2, 2.0 * rh, sh2)

    return S, Sc


def simulate_quantised(S0, r, sigma, T, M, n_paths, rng, F):
    h = T / M
    rh = F.rh.const(r * h)
    sh2 = F.sh2.const(sigma * np.sqrt(h))
    S = F.S.const(np.full(n_paths, S0))
    for _ in range(M):
        S = step_fixed(S, F.dW.const(rng.standard_normal(n_paths)), rh, sh2, F)
    return S


def payoff(S):
    return np.maximum(S - 1.0, 0.0)


def _price(terminal, K, r, T):
    p = payoff(terminal)
    scale = K * np.exp(-r * T)
    return scale * p.mean(), scale * p.std(ddof=1) / np.sqrt(p.size)


def mc_price(S0, K, r, sigma, T, M, n_paths, rng):
    return _price(simulate(S0 / K, r, sigma, T, M, n_paths, rng), K, r, T)


def mc_price_quantised(S0, K, r, sigma, T, M, n_paths, rng, F):
    return _price(simulate_quantised(S0 / K, r, sigma, T, M, n_paths, rng, F).to_float(), K, r, T)


def level_delta(S0, K, r, sigma, T, level, n_paths, rng):
    scale = K * np.exp(-r * T)
    if level == 0:
        return scale * payoff(simulate(S0 / K, r, sigma, T, 1, n_paths, rng))
    S, Sc = simulate_coupled(S0 / K, r, sigma, T, 2 ** level, n_paths, rng)
    return scale * (payoff(S) - payoff(Sc))
