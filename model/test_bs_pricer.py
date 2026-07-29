import numpy as np

import bs_pricer_model as bs
from checks import Checks

check = Checks()

c = bs.call(100.0, 100.0, 0.05, 0.2, 1.0)
p = bs.put(100.0, 100.0, 0.05, 0.2, 1.0)
check("textbook ATM call", abs(c - 10.450583572185565) < 1e-12, f"got {c}")
check("textbook ATM put", abs(p - 5.573526022256971) < 1e-12, f"got {p}")

rng = np.random.default_rng(20260729)
n = 20000
S = rng.uniform(10.0, 300.0, n)
K = rng.uniform(10.0, 300.0, n)
r = rng.uniform(-0.02, 0.10, n)
sg = rng.uniform(0.01, 1.50, n)
T = rng.uniform(1e-3, 5.0, n)

res = np.abs(bs.parity_residual(S, K, r, sg, T))
check("put-call parity", res.max() < 1e-9, f"max residual {res.max():.3e}")

scaled = np.abs(bs.call(S, K, r, sg, T) - K * bs.call_normalised(S / K, r, sg, T))
check("K-homogeneity", scaled.max() < 1e-9, f"max dev {scaled.max():.3e}")

c_all = bs.call(S, K, r, sg, T)
lower = np.maximum(S - K * np.exp(-r * T), 0.0)
check("no-arbitrage bounds", bool(np.all(c_all >= lower - 1e-9) and np.all(c_all <= S + 1e-9)))

# Both degenerate paths must land on discounted intrinsic, not NaN.
z_sig = bs.call(np.array([120.0, 80.0]), 100.0, 0.05, 0.0, 1.0)
z_exp = bs.call(np.array([120.0, 80.0]), 100.0, 0.05, 0.2, 0.0)
check("sigma -> 0 degenerate", np.allclose(z_sig, [120.0 - 100.0 * np.exp(-0.05), 0.0]), f"got {z_sig}")
check("T -> 0 degenerate", np.allclose(z_exp, [20.0, 0.0]), f"got {z_exp}")

g = bs.greeks(100.0, 100.0, 0.05, 0.2, 1.0)
h = 1e-5
fd = (bs.call(100.0 + h, 100.0, 0.05, 0.2, 1.0) - bs.call(100.0 - h, 100.0, 0.05, 0.2, 1.0)) / (2 * h)
check("delta vs finite difference", abs(g["delta"] - fd) < 1e-7, f"{g['delta']:.9f} vs {fd:.9f}")

check.report()
