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
    a = pad_coeffs(coeffs)

    x2 = x * x
    x4 = x2 * x2

    t0 = a[0] + a[1] * x
    t1 = a[2] + a[3] * x
    t2 = a[4] + a[5] * x

    return (t0 + x2 * t1) + x4 * t2


def evaluate_estrin_fixed(coeffs, x, coef_fmt, x_fmt, work_fmt):
    a = []
    for c in pad_coeffs(coeffs):
        a.append(coef_fmt.const(c))

    xf = x_fmt.const(x)

    x2 = work_fmt.cast(xf * xf)
    x4 = work_fmt.cast(x2 * x2)

    t0 = work_fmt.cast(a[0] + work_fmt.cast(a[1] * xf))
    t1 = work_fmt.cast(a[2] + work_fmt.cast(a[3] * xf))
    t2 = work_fmt.cast(a[4] + work_fmt.cast(a[5] * xf))

    left = work_fmt.cast(t0 + work_fmt.cast(x2 * t1))
    right = work_fmt.cast(x4 * t2)
    return work_fmt.cast(left + right)


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
