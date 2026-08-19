import numpy as np


def _rne(code, shift):
    # Round to nearest even. Unbiased, so it is what a serial recurrence wants,
    # but it costs a comparator and an LSB test in hardware.
    if shift <= 0:
        return code << -shift
    q = code >> shift
    rem = code - (q << shift)
    half = 1 << (shift - 1)
    return q + ((rem > half) | ((rem == half) & (q & 1).astype(bool)))


def _half_up(code, shift):
    # Round half away from zero upward: one adder. Biased by half an LSB, which
    # only matters where roundings accumulate in the same direction.
    if shift <= 0:
        return code << -shift
    return (code + (1 << (shift - 1))) >> shift


ROUNDING = {"rne": _rne, "half_up": _half_up}


class Fx:
    # A value on a wire: an integer code and the scale it is understood at.
    def __init__(self, code, frac):
        self.code, self.frac = np.asarray(code, dtype=np.int64), frac

    def __mul__(self, other):
        return Fx(self.code * other.code, self.frac + other.frac)

    def __add__(self, other):
        f = max(self.frac, other.frac)
        return Fx((self.code << (f - self.frac)) + (other.code << (f - other.frac)), f)

    def to_float(self):
        return np.asarray(self.code, dtype=float) * (2.0 ** -self.frac)


class Fmt:
    # register: how many bits, where its binary point sits, how it rounds.
    def __init__(self, width, frac, rounding="rne"):
        if not 1 <= width <= 31:
            # Products are held in int64, so two operands must fit in 62 bits.
            raise ValueError(f"width {width} outside 1..31")
        if rounding not in ROUNDING:
            raise ValueError(f"unknown rounding {rounding}")
        self.width, self.frac, self.rounding = width, frac, rounding
        self.lo = -(1 << (width - 1))
        self.hi = (1 << (width - 1)) - 1
        self.saturations = 0

    @classmethod
    def for_range(cls, width, max_abs, rounding="rne"):
        # floor(log2) is the position of the leading bit, +1 for that bit and
        # +1 for the sign. ceil() is wrong at exact powers of two: it gives 3
        # integer bits for 4.0, whose largest value is 3.99988, so 4.0 saturates.
        int_bits = 1 if max_abs <= 0 else int(np.floor(np.log2(max_abs))) + 2
        return cls(width, width - int_bits, rounding)

    def _clip(self, code):
        self.saturations += int(np.count_nonzero((code < self.lo) | (code > self.hi)))
        return np.clip(code, self.lo, self.hi)

    def cast(self, x):
        return Fx(self._clip(ROUNDING[self.rounding](x.code, x.frac - self.frac)), self.frac)

    def const(self, value):
        scaled = np.asarray(value, dtype=float) * (2.0 ** self.frac)
        return Fx(self._clip(np.rint(scaled).astype(np.int64)), self.frac)
