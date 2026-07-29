import numpy as np

import bs_pricer_model as bs
from checks import Checks
from fixedpoint import Formats
from path_engine_model import level_delta, mc_price, mc_price_quantised
from path_step_model import measure_ranges


def slope(x, y):
    return np.polyfit(np.log(x), np.log(y), 1)[0]


check = Checks()

S0, K, r, sigma, T = 100.0, 100.0, 0.05, 0.2, 1.0
exact = bs.call(S0, K, r, sigma, T)
rng = np.random.default_rng(20260729)
print(f"exact BS call = {exact:.6f}\n")

px, se = mc_price(S0, K, r, sigma, T, 256, 2_000_000, rng)
check("MC brackets exact price", abs(px - exact) < 3 * se + 0.01,
      f"{px:.4f} +/- {se:.4f}  (err {px - exact:+.4f})")

counts = np.array([2000, 8000, 32000, 128000, 512000])
spread = [np.std([mc_price(S0, K, r, sigma, T, 16, int(n), rng)[0] for _ in range(12)], ddof=1)
          for n in counts]
s = slope(counts, np.array(spread))
check("statistical slope -0.5", -0.65 < s < -0.35, f"fitted {s:+.3f}")

corr = np.array([abs(level_delta(S0, K, r, sigma, T, l, 2_000_000, rng).mean())
                 for l in (3, 4, 5)])
ratios = corr[:-1] / corr[1:]
check("weak convergence alpha ~ 1", bool(np.all((ratios > 1.8) & (ratios < 2.6))),
      "ratios " + " ".join(f"{x:.2f}" for x in ratios))

lv = np.array([1, 2, 3, 4, 5])
var = [level_delta(S0, K, r, sigma, T, int(l), 200_000, rng).var(ddof=1) for l in lv]
check("level variance beta ~ 1", 0.7 < slope(T / 2.0 ** lv, np.array(var)) < 1.35,
      f"fitted {slope(T / 2.0 ** lv, np.array(var)):+.3f}")

tel = sum(level_delta(S0, K, r, sigma, T, l, 400_000, rng).mean() for l in range(6))
check("telescoping sum reaches price", abs(tel - exact) < 0.05, f"{tel:.4f} vs {exact:.4f}")

M = 32
ranges = measure_ranges(1.0, r, sigma, T, M, 200_000, np.random.default_rng(5))
ref, _ = mc_price(S0, K, r, sigma, T, M, 400_000, np.random.default_rng(99))
print(f"\n{'width':>6} {'price':>10} {'vs float':>10}")
errs = []
for w in (8, 10, 12, 16, 20, 24):
    p, _ = mc_price_quantised(S0, K, r, sigma, T, M, 400_000,
                              np.random.default_rng(99), Formats(w, ranges))
    errs.append(abs(p - ref))
    print(f"{w:6d} {p:10.5f} {errs[-1]:10.5f}")
print()

check("narrow widths cost accuracy", errs[0] > 10 * errs[-1],
      f"8-bit {errs[0]:.4f} vs 24-bit {errs[-1]:.5f}")
check("wide widths track the float engine", errs[-1] < 0.01, f"24-bit err {errs[-1]:.5f}")
check("error shrinks monotonically in width", all(np.diff(errs) < 1e-9))

check.report()
