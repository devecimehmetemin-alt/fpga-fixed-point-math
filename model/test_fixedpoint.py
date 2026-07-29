import numpy as np

from checks import Checks
from fixedpoint import Fmt, Fx

check = Checks()

q8 = Fmt(8, 0)
ties = Fx(np.array([1, 3, 5, -1, -3, 7]), 1)
check("round-to-nearest-even", q8.cast(ties).code.tolist() == [0, 2, 2, 0, -2, 4],
      f"0.5,1.5,2.5,-0.5,-1.5,3.5 -> {q8.cast(ties).code.tolist()}")

a, b = Fmt(16, 8).const(1.5), Fmt(16, 12).const(0.0625)
check("addition is exact", (a + b).to_float() == 1.5625, f"got {(a + b).to_float()}")

p = Fmt(16, 8).const(1.5) * Fmt(16, 8).const(2.5)
check("multiply is exact", p.to_float() == 3.75, f"got {p.to_float()}")

sat = Fmt(8, 4)
sat.const(np.array([7.9, 100.0, -100.0]))
check("saturation is counted", sat.saturations == 2, f"counted {sat.saturations}")

f = Fmt.for_range(16, 3.07)
check("range picks the binary point", (f.width, f.frac) == (16, 13), f"got w={f.width} f={f.frac}")

q = Fmt(12, 8)
x = np.random.default_rng(1).uniform(-7.0, 7.0, 10000)
check("round trip within half an LSB", np.abs(q.const(x).to_float() - x).max() <= 2.0 ** -9 + 1e-15)

check.report()
