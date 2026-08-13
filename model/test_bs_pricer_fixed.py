"""The fixed-point pricer against the float reference, over the declared box."""

import math

import numpy as np

import bs_pricer_fixed as fx
import bs_pricer_model as bs
from checks import Checks

check = Checks()
fmt = fx.DEFAULT
rng = np.random.default_rng(20260811)
n = 4000

S, K, r, sg, T = bs.sample(n, rng)
rows = [fx.encode(*t, fmt=fmt) for t in zip(S, K, r, sg, T)]
got = fx.evaluate_float(rows, fmt)
gc = np.array([c for c, _ in got])
gp = np.array([p for _, p in got])

# The reference is the clamped price: the clamp is the spec, not an error.
ref_c = bs.call(S, K, r, sg, T, clamp=bs.D_CLAMP)
ref_p = bs.put(S, K, r, sg, T, clamp=bs.D_CLAMP)

# Relative to K, because the datapath is K-normalised and every absolute error
# in the interior is scaled by K exactly once at the end.
err_c = np.abs(gc - ref_c) / K
err_p = np.abs(gp - ref_p) / K
bits = lambda v: float("inf") if v <= 0 else -math.log2(v)

v = sg * np.sqrt(T)
check("call error under 2^-9 of K worst case", err_c.max() < 2.0 ** -9,
      f"max {err_c.max():.3e} ({bits(err_c.max()):.1f} bits) at v={v[err_c.argmax()]:.2e}")
check("put error under 2^-9 of K worst case", err_p.max() < 2.0 ** -9,
      f"max {err_p.max():.3e} ({bits(err_p.max()):.1f} bits)")

# The typical case is what the module is actually for. v > 0.05 is sigma and T
# both away from their floors.
mid = v > 0.05
check("call error under 2^-16 of K for v > 0.05", err_c[mid].max() < 2.0 ** -16,
      f"max {err_c[mid].max():.3e} ({bits(err_c[mid].max()):.1f} bits), "
      f"{mid.sum()}/{n} samples")

# Parity is a statement the pricer makes about itself, so it holds even where
# the price is inaccurate. This is the check that survives onto the RTL with no
# reference model at all.
res = np.abs(fx.parity_residual(rows, fmt)) / K
check("put-call parity on the module's own numbers", res.max() < 2.0 ** -18,
      f"max {res.max():.3e} ({bits(res.max()):.1f} bits)")

# No output may be negative: a call is worth at least nothing.
check("no negative prices", bool(np.all(gc >= 0) and np.all(gp >= 0)),
      f"min call {gc.min():.6g}, min put {gp.min():.6g}")

# Both saturating paths must be reachable, or the corner tests are decorative.
s = [fx.stages(*row, fmt=fmt) for row in rows]
sat = sum(1 for e in s if abs(e["d1_wide"]) > fmt.CLAMP_CODE)
check("clamp is exercised", sat > 0, f"{sat}/{n} samples clamp d1")

# The error is set by x = S/K, not by 1/v. Both claims get pinned, because the
# conditioning story is the one everybody reaches for first and it is wrong
# here: wherever 1/v is large, d is past the clamp and N has already saturated.
xf = np.array([s["x"] for s in [fx.stages(*r, fmt=fmt) for r in rows]])
xf = xf * 2.0 ** -fmt.X_F
check("error grows with x = S/K",
      err_c[xf > 10].max() > 2 * err_c[xf < 1].max(),
      f"x>10 max {err_c[xf > 10].max():.2e}, x<1 max {err_c[xf < 1].max():.2e}, "
      f"ratio {err_c[xf > 10].max() / err_c[xf < 1].max():.1f}x")
check("error does NOT grow as v -> 0",
      err_c[v < 0.01].max() < err_c[v > 0.5].max(),
      f"v<0.01 max {err_c[v < 0.01].max():.2e}, "
      f"v>0.5 max {err_c[v > 0.5].max():.2e}")

check.report()

print()
print(f"{'v = sigma*sqrt(T)':<22} {'samples':>8} {'max err / K':>13} {'bits':>7}")
edges = [0.0, 0.003, 0.01, 0.03, 0.1, 0.3, 1.0, 4.0]
for a, b in zip(edges[:-1], edges[1:]):
    sel = (v >= a) & (v < b)
    if sel.sum():
        e = err_c[sel].max()
        print(f"  [{a:6.3f}, {b:6.3f})       {sel.sum():>8} {e:>13.3e} {bits(e):>7.1f}")
