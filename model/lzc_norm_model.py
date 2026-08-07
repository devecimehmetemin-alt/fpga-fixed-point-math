"""Golden model for rtl/lzc_norm.sv.

Splits an unsigned integer into m * 2^e with m in [1, 2)
the fixed-to-floating decomposition the ln/sqrt/recip range reductions need.

e = floor(log2 x) is the position of the highest set bit, so a leading-zero
count answers it exactly. Nothing here approximates: x == m * 2^e holds with
equality for every input.

x is a raw W-bit unsigned integer. A caller holding Q(W-F).F must use e - F as
the true exponent; this model deliberately does not know F.
"""


def normalize(x, W):
    """Returns (m, e, lz, zero) -- the outputs the RTL must produce."""
    assert 0 <= x < (1 << W), f"x={x} does not fit in {W} bits"
    if x == 0:
        # No m >= 1 gives m * 2^e == 0. Pin the outputs so every input has a
        # checkable value, and flag it. lz = W is the true count, but the RTL
        # saturates one short -- do not compare lz here.
        return 0, 0, W, True
    lz = W - x.bit_length()
    return x << lz, W - 1 - lz, lz, False


def corners(W):
    """Every exponent boundary: 2^k is the first x with lz = W-1-k, 2^k - 1
    the last with one more leading zero. Straddles every layer transition."""
    xs = {0, (1 << W) - 1}
    for k in range(W):
        xs |= {(1 << k) - 1, 1 << k, (1 << k) + 1}
    return sorted(v for v in xs if 0 <= v < (1 << W))


if __name__ == "__main__":
    W = 12
    for x in range(1 << W):
        m, e, lz, zero = normalize(x, W)
        assert (m >> lz if not zero else 0) == x  # x == m * 2^e, exactly
        assert zero or m >> (W - 1) == 1  # m is normalised
    print(f"W={W}: x == m*2^e for all {1 << W} inputs")

    m, e, lz, zero = normalize(22, 8)
    print(f"x=22 W=8 -> m={m:08b} = {m / 128}, e={e}, lz={lz}")
