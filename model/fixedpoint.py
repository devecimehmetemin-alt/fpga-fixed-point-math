import numpy as np

VARS = ("dW", "sh2", "rh", "mult1", "sum1", "mult2", "S", "mult2c", "Sc")


def _rne(code, shift):
    # Round to nearest even
    if shift <= 0:
        return code << -shift
    q = code >> shift
    rem = code - (q << shift)
    half = 1 << (shift - 1)
    return q + ((rem > half) | ((rem == half) & (q & 1).astype(bool)))


class Fx:
    # A value on a wire: an integer code and the scale it is understood at.
    __slots__ = ("code", "frac")

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
    # register: how many bits, and where its binary point sits.
    def __init__(self, width, frac, name=""):
        if not 1 <= width <= 31:
            # Products are held in int64, so two operands must fit in 62 bits.
            raise ValueError(f"width {width} outside 1..31")
        self.width, self.frac, self.name = width, frac, name
        self.lo = -(1 << (width - 1))
        self.hi = (1 << (width - 1)) - 1
        self.saturations = 0

    @classmethod
    def for_range(cls, width, max_abs, name=""):
        int_bits = 1 if max_abs <= 0 else int(np.ceil(np.log2(max_abs))) + 1
        return cls(width, width - int_bits, name)

    def _clip(self, code):
        self.saturations += int(np.count_nonzero((code < self.lo) | (code > self.hi)))
        return np.clip(code, self.lo, self.hi)

    def cast(self, x):
        return Fx(self._clip(_rne(x.code, x.frac - self.frac)), self.frac)

    def const(self, value):
        scaled = np.asarray(value, dtype=float) * (2.0 ** self.frac)
        return Fx(self._clip(np.rint(scaled).astype(np.int64)), self.frac)

    def __repr__(self):
        return f"Fmt({self.name or '?'}, w={self.width}, f={self.frac})"


class Formats:
    def __init__(self, widths, ranges):
        if isinstance(widths, int):
            widths = {v: widths for v in VARS}
        for v in VARS:
            setattr(self, v, Fmt.for_range(widths[v], ranges[v], v))

    def saturations(self):
        return {v: getattr(self, v).saturations for v in VARS if getattr(self, v).saturations}
