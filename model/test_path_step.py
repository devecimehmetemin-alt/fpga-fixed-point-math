import numpy as np

from checks import Checks
from fixedpoint import Formats
from path_step_model import measure_ranges, step_fixed, step_float

check = Checks()

r, sigma, T, M = 0.05, 0.2, 1.0, 32
h = T / M
rh_f, sh2_f = r * h, sigma * np.sqrt(h)

ranges = measure_ranges(1.0, r, sigma, T, M, 100_000, np.random.default_rng(5))
print("ranges: " + "  ".join(f"{k}={v:.3g}" for k, v in ranges.items()) + "\n")

z = np.random.default_rng(1).standard_normal(50_000)
S0 = np.full(z.size, 1.0)
ref = step_float(S0, z, rh_f, sh2_f)

print(f"{'width':>6} {'max |err|':>12}")
errs = []
for w in (8, 12, 16, 20, 24):
    F = Formats(w, ranges)
    got = step_fixed(F.S.const(S0), F.dW.const(z), F.rh.const(rh_f),
                     F.sh2.const(sh2_f), F).to_float()
    errs.append(float(np.abs(got - ref).max()))
    print(f"{w:6d} {errs[-1]:12.3e}")
print()

check("one step tracks the float model", errs[-1] < 1e-6, f"24-bit max err {errs[-1]:.2e}")
check("error falls monotonically with width", all(np.diff(errs) < 0))
check("roughly one bit per halving", 8 < errs[1] / errs[2] < 32, f"12->16 bits: {errs[1] / errs[2]:.1f}x")

F = Formats(16, ranges)
args = (F.S.const(S0), F.dW.const(z), F.rh.const(rh_f), F.sh2.const(sh2_f), F)
check("bit-exact reproducibility", np.array_equal(step_fixed(*args).code, step_fixed(*args).code))
check("no silent saturation at w=16", not F.saturations(), f"{F.saturations()}")

check.report()
