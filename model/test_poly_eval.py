import numpy as np

import poly_eval_model as poly
from checks import Checks
from fixedpoint import Fmt

check = Checks()

coeffs = [0.5, -1.25, 2.0, 0.75, -0.125, 0.0625]
x = np.linspace(-0.5, 0.5, 1001)

h = poly.evaluate_horner(coeffs, x)
e = poly.evaluate_estrin(coeffs, x)
check("Estrin matches Horner", np.abs(h - e).max() < 1e-12, f"max diff {np.abs(h - e).max():.2e}")

exp2 = poly.exp2_coefficients()
ln1p = poly.ln1p_coefficients()
err_exp2 = poly.max_error(lambda f: 2.0 ** f, exp2, -0.5, 0.5)
err_ln1p = poly.max_error(lambda z: np.log(1.0 + z), ln1p, 0.0, 1.0)

print("2^f  coefficients: " + " ".join(f"{c:+.6f}" for c in exp2))
print("ln1p coefficients: " + " ".join(f"{c:+.6f}" for c in ln1p))
print(f"\n2^f  max error over [-0.5, 0.5] = {err_exp2:.3e}")
print(f"ln1p max error over [0, 1)     = {err_ln1p:.3e}\n")

check("degree-5 fit of 2^f reaches single precision", err_exp2 < 1e-7, f"{err_exp2:.2e}")
check("degree-5 fit of ln(1+z) is usable", err_ln1p < 1e-4, f"{err_ln1p:.2e}")

check("2^f is the easier fit", err_exp2 < err_ln1p)

ref = poly.evaluate_estrin(exp2, x)
print(f"{'width':>6} {'max |err|':>12}")
errs = []
for w in (10, 12, 16, 20, 24):
    coef_fmt = Fmt(w, w - 2)
    x_fmt = Fmt(w, w - 2)
    work_fmt = Fmt(w, w - 3)
    got = poly.evaluate_estrin_fixed(exp2, x, coef_fmt, x_fmt, work_fmt).to_float()
    errs.append(float(np.abs(got - ref).max()))
    print(f"{w:6d} {errs[-1]:12.3e}")
print()

check("fixed point tracks float as width grows", errs[-1] < 1e-5, f"24-bit {errs[-1]:.2e}")
check("error falls monotonically with width", all(np.diff(errs) < 0))

check.report()
