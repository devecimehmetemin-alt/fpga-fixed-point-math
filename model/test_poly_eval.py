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

# Datapath sweep. The reference is the same polynomial in float, so this is
# quantisation error alone -- the fit error above is excluded. Feeds the three
# open decisions on poly_eval: working width, rounding mode, and whether the
# saturation logic is needed at all.
cases = [
    ("2^f", exp2, -0.5, 0.5),
    ("ln1p", ln1p, 0.0, 1.0),
]
widths = (16, 18, 20, 24)
modes = ("rne", "half_up")

print(f"{'fit':>5} {'round':>8} {'W_W':>4} {'max |err|':>12} {'sat':>5}")
results = {}
for name, cf, lo, hi in cases:
    xs = np.linspace(lo, hi, 4001)
    ref = poly.evaluate_estrin(cf, xs)
    for mode in modes:
        for w in widths:
            coef_fmt = Fmt(w, w - 2, mode)
            x_fmt = Fmt(w, w - 2, mode)
            work_fmt = Fmt(w, w - 3, mode)
            got = poly.evaluate_estrin_fixed(cf, xs, coef_fmt, x_fmt, work_fmt).to_float()
            err = float(np.abs(got - ref).max())
            sat = coef_fmt.saturations + x_fmt.saturations + work_fmt.saturations
            results[(name, mode, w)] = (err, sat)
            print(f"{name:>5} {mode:>8} {w:4d} {err:12.3e} {sat:5d}")
print()

for name, _, _, _ in cases:
    for mode in modes:
        errs = [results[(name, mode, w)][0] for w in widths]
        check(f"{name:>4} {mode:>8}: error falls with width", all(np.diff(errs) < 0))

sats = sum(results[k][1] for k in results)
check("no saturation anywhere in the sweep", sats == 0, f"{sats} clipped values")

for name, _, _, _ in cases:
    rne = results[(name, "rne", 18)][0]
    hup = results[(name, "half_up", 18)][0]
    check(f"{name:>4}: half_up within 2x of rne at 18 bits", hup < 2.0 * rne,
          f"rne {rne:.2e}  half_up {hup:.2e}")

check.report()
