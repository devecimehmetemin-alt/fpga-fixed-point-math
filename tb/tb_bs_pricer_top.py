"""cocotb testbench for bs_pricer_top. Drive the box corners plus a random
sample, compare both prices to the golden model bit for bit, then check the
properties a bit-exact compare cannot see: non-negativity, put-call parity, and
monotonicity in the spot."""

import itertools
import math
import os

import cocotb
import numpy as np
from cocotb.clock import Clock
from cocotb.triggers import ClockCycles, FallingEdge, RisingEdge

import bs_pricer_fixed as fx
import bs_pricer_model as bs

SK_W, SK_F = int(os.environ["SK_W"]), int(os.environ["SK_F"])
T_W, T_F = int(os.environ["T_W"]), int(os.environ["T_F"])
SG_W, SG_F = int(os.environ["SG_W"]), int(os.environ["SG_F"])
R_W, R_F = int(os.environ["R_W"]), int(os.environ["R_F"])
OUT_F = int(os.environ["OUT_F"])
N_F = int(os.environ["N_F"])
LATENCY = int(os.environ["LATENCY"])

NAMES = ("SK_W", "SK_F", "T_W", "T_F", "SG_W", "SG_F", "R_W", "R_F",
         "IK_W", "IK_F", "X_W", "X_F", "LX_W", "LX_F", "ST_W", "ST_F",
         "V_W", "V_F", "IV_W", "IV_F", "HV_W", "HV_F", "HT_W", "HT_F",
         "NM_W", "NM_F", "D_W", "IN_W", "IN_F", "U_W", "U_F", "DC_W", "DC_F",
         "CN_W", "CN_F", "PN_W", "PN_F", "OUT_W", "OUT_F",
         "K_IDX", "V_IDX", "EC_W", "EC_F", "EW_W", "EW_F")

FMT = fx.PricerFormats(**{n: int(os.environ[n]) for n in NAMES})
LSB = 2.0 ** -OUT_F


def mask(value, width):
    # A Python int to wire
    return int(value) & ((1 << width) - 1)


async def run(dut, rows):
    """Drive rows of input codes back to back, return one (call, put) per row."""
    cocotb.start_soon(Clock(dut.clk, 10, "ns").start())
    dut.rst.value = 1
    dut.in_valid.value = 0
    for sig in (dut.s, dut.k, dut.r, dut.sigma, dut.tau):
        sig.value = 0

    await ClockCycles(dut.clk, 3)
    dut.rst.value = 0

    got = []
    for i in range(len(rows) + LATENCY + 2):
        await RisingEdge(dut.clk)
        live = i < len(rows)
        s, k, r, sg, t = rows[i] if live else (0, 0, 0, 0, 0)
        dut.in_valid.value = int(live)
        dut.s.value = mask(s, SK_W)
        dut.k.value = mask(k, SK_W)
        dut.r.value = mask(r, R_W)
        dut.sigma.value = mask(sg, SG_W)
        dut.tau.value = mask(t, T_W)

        await FallingEdge(dut.clk)
        if int(dut.out_valid.value):
            got.append((dut.call.value.to_signed(), dut.put.value.to_signed()))
    return got


def stimulus(n, seed):
    """Box corners first, then a random fill. The corners are the samples the
    uniform draw will never produce and the saturating paths live on them."""
    corners = list(itertools.product(*[bs.RANGES[a] for a in bs.ARGS]))
    rng = np.random.default_rng(seed)
    drawn = list(zip(*bs.sample(n, rng)))
    return corners + [tuple(float(v) for v in row) for row in drawn]


@cocotb.test()
async def test_box(dut):
    reals = stimulus(2000, 20260813)
    rows = [fx.encode(*args, fmt=FMT) for args in reals]

    golden = fx.evaluate_fixed(rows, FMT)
    got = await run(dut, rows)

    assert len(got) == len(golden), \
        f"got {len(got)} results, expected {len(golden)} -- LATENCY is wrong"

    for i, ((gc, gp), (ec, ep)) in enumerate(zip(got, golden)):
        S, K, r, sg, T = reals[i]
        assert (gc, gp) == (ec, ep), (
            f"sample {i}: S={S:.6f} K={K:.6f} r={r:+.6f} sigma={sg:.6f} T={T:.6f}\n"
            f"  call rtl={gc} model={ec} diff={gc - ec}\n"
            f"  put  rtl={gp} model={ep} diff={gp - ep}")

    # Structural, on the RTL's own outputs. An option is never worth less than
    # nothing, and the floor at the last stage is the only thing enforcing it.
    neg = [(i, c, p) for i, (c, p) in enumerate(got) if c < 0 or p < 0]
    assert not neg, f"{len(neg)} negative prices, first at sample {neg[0]}"

    # Put-call parity: C - P = S - K*e^(-rT). This is the check that survives
    # onto silicon with no model at all. Samples where either price floored are
    # excluded -- the floor is exactly where parity stops holding, by design.
    # Together with C, P >= 0 this also pins the no-arbitrage bounds, so they
    # are not checked separately.
    res, live = [], 0
    for (gc, gp), (S, K, r, sg, T) in zip(got, reals):
        if gc == 0 or gp == 0:
            continue
        live += 1
        res.append(abs((gc - gp) * LSB - (S - K * math.exp(-r * T))) / K)
    worst = max(res)
    assert worst < 2.0 ** -18, (
        f"parity residual {worst:.3e} of K ({-math.log2(worst):.1f} bits) "
        f"over {live} unfloored samples")

    dut._log.info(
        f"{len(got)} samples bit-exact, 0 negative prices, "
        f"parity {-math.log2(worst):.1f} bits over {live} unfloored samples")


@cocotb.test()
async def test_monotone_in_spot(dut):
    """C is non-decreasing in S and P is non-increasing, at everything else
    fixed. Invisible to a bit-exact compare, which would pass just as happily if
    the model stepped backwards too, and it is the property a caller doing an
    implied-vol solve on top of this actually leans on."""
    K, r, sg, T = 100.0, 0.03, 0.25, 1.0
    lo, hi = bs.RANGES["S"]
    count = 801
    reals = [(lo + (hi - lo) * i / (count - 1), K, r, sg, T)
             for i in range(count)]
    rows = [fx.encode(*args, fmt=FMT) for args in reals]

    got = await run(dut, rows)
    assert len(got) == count

    calls = [c for c, _ in got]
    puts = [p for _, p in got]

    # Both prices are K*(a*N(d1) -+ b*N(d2)) with a = S/K and b = disc, so one
    # code of either CDF moves the output by (S + K*disc)*2^(OUT_F-N_F) lsb.
    # Deep out of the money the two terms are nearly equal and each steps by
    # more than their difference, so the price is flat to the eye and jitters by
    # a code. That is N's output format showing through, not an inversion, and
    # the only lever on it is N_F. Above the grain the ordering has to be exact.
    def grain(S):
        return int(math.ceil((S + K * math.exp(-r * T)) * 2.0 ** (OUT_F - N_F))) + 2

    for name, seq, sign in (("call", calls, +1), ("put", puts, -1)):
        bad = [(reals[i][0], seq[i], seq[i + 1], grain(reals[i + 1][0]))
               for i in range(count - 1) if sign * (seq[i + 1] - seq[i]) < 0]
        past = [e for e in bad if abs(e[2] - e[1]) > e[3]]
        assert not past, (
            f"{name} reverses past what one code of N is worth, {len(past)} "
            f"times, first at S={past[0][0]:.4f}: {past[0][1]} -> {past[0][2]} "
            f"against a grain of {past[0][3]}")

        # Where the price is well clear of the grain the two terms are no longer
        # cancelling and the ordering is a hard property, not a statistical one.
        strict = [e for e in bad if max(e[1], e[2]) > 16 * e[3]]
        assert not strict, (
            f"{name} reverses at {strict[0][1]} lsb, far above its "
            f"{strict[0][3]} lsb grain, {len(strict)} times")
        worst = max((abs(b - a) for _, a, b, _ in bad), default=0)
        dut._log.info(
            f"{name}: {len(bad)} reversals, worst {worst} lsb, "
            f"grain {grain(lo)}-{grain(hi)} lsb across the sweep")

    # The slope of C in S is N(d1), which is in (0,1), so the price may sit flat
    # over an lsb but must not outrun its own argument.
    span = (calls[-1] - calls[0]) * LSB
    assert span <= (hi - lo) + LSB,         f"call rose {span:.6f} over an S range of {hi - lo:.6f}"

    dut._log.info(f"{count} spots, call span {span:.4f} of {hi - lo:.4f}")
