import numpy as np

DEGREE = 5
N_COEFFS = DEGREE + 1


def pad_coeffs(coeffs):
    a = list(coeffs)
    while len(a) < N_COEFFS:
        a.append(0.0)
    return a


def evaluate_horner(coeffs, x):
    result = 0.0
    for c in reversed(coeffs):
        result = result * x + c
    return result


def evaluate_estrin(coeffs, x):
    # p = t0 + x2*(t1 + x2*t2). Six multiplies, arithmetic depth 3.
    a = pad_coeffs(coeffs)

    x2 = x * x

    t0 = a[0] + a[1] * x
    t1 = a[2] + a[3] * x
    t2 = a[4] + a[5] * x

    return t0 + x2 * (t1 + x2 * t2)


def evaluate_estrin_fixed(coeffs, x, coef_fmt, x_fmt, work_fmt):
    # Mirrors rtl/poly_eval.sv one for one. Six multiplies, and exactly one
    # cast per multiply because a DSP48 carries the full product into its own
    # post-adder and only the sum leaves the block. Rounding anywhere else
    # would model a circuit that isn't there.
    a = []
    for c in pad_coeffs(coeffs):
        a.append(coef_fmt.const(c))

    xf = x_fmt.const(x)

    # pipeline stages 1-2: four products, three multiply-adds
    x2 = work_fmt.cast(xf * xf)
    t0 = work_fmt.cast(a[0] + a[1] * xf)
    t1 = work_fmt.cast(a[2] + a[3] * xf)
    t2 = work_fmt.cast(a[4] + a[5] * xf)

    # stages 3-4
    w = work_fmt.cast(t1 + x2 * t2)

    # stages 5-6
    return work_fmt.cast(t0 + x2 * w)


def chebyshev_nodes(lo, hi, count):
    k = np.arange(count)
    unit = np.cos((2.0 * k + 1.0) * np.pi / (2.0 * count))
    return 0.5 * (lo + hi) + 0.5 * (hi - lo) * unit


def fit_coefficients(func, lo, hi, degree=DEGREE, count=200):
    x = chebyshev_nodes(lo, hi, count)
    y = func(x)
    fitted = np.polyfit(x, y, degree)
    # polyfit gives highest power first, the RTL wants a0 first
    return fitted[::-1]


def max_error(func, coeffs, lo, hi, count=20000):
    x = np.linspace(lo, hi, count)
    return float(np.abs(evaluate_horner(coeffs, x) - func(x)).max())


def exp2_coefficients():
    # 2^f on [-0.5, 0.5], for exp_unit after range reduction
    return fit_coefficients(lambda f: 2.0 ** f, -0.5, 0.5)


def ln1p_coefficients():
    # ln(1+z) on [0, 1), for ln_unit after the mantissa split
    return fit_coefficients(lambda z: np.log(1.0 + z), 0.0, 1.0)
