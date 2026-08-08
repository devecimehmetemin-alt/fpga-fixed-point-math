"""cocotb testbench for ln_unit. drive a sweep, compare every output to the
golden model bit for bit."""

import os

import numpy as np

import cocotb
from cocotb.clock import Clock
from cocotb.triggers import ClockCycles, FallingEdge, RisingEdge

import ln_unit_model as ln
import lzc_norm_model as lzc

W, F = int(os.environ["W"]), int(os.environ["F"])
X_W, X_F = int(os.environ["X_W"]), int(os.environ["X_F"])
C_W, C_F = int(os.environ["C_W"]), int(os.environ["C_F"])
W_W, W_F = int(os.environ["W_W"]), int(os.environ["W_F"])
L_W, Y_W = int(os.environ["L_W"]), int(os.environ["Y_W"])

LATENCY = int(os.environ["LATENCY"])

FMT = ln.LnFormats(W=W, F=F, X_W=X_W, X_F=X_F, C_W=C_W, C_F=C_F,
                   W_W=W_W, W_F=W_F, L_W=L_W, Y_W=Y_W)


def mask(value, width):
    # A Python int to wire
    return int(value) & ((1 << width) - 1)


def stimulus():
    # corners is testing for transitions. 2^k is the first x with a
    # given e, 2^k - 1 the last with one less.
    xs = lzc.corners(W)
    # S/K near 1, where e - F crosses from -1 to 0 and prod changes sign.
    xs += [int(round(r * 2 ** F)) for r in np.linspace(0.5, 2.0, 257)]
    xs += [int(v) for v in np.random.default_rng(0).integers(0, 1 << W, 512)]
    return xs


@cocotb.test()
async def test_sweep(dut):
    xs = stimulus()
    golden = [(int(y), int(x == 0)) for y, x in zip(ln.evaluate_fixed(xs, FMT), xs)]

    cocotb.start_soon(Clock(dut.clk, 10, "ns").start())
    dut.rst.value = 1
    dut.in_valid.value = 0
    dut.x.value = 0

    await ClockCycles(dut.clk, 3)
    dut.rst.value = 0

    got = []
    for i in range(len(xs) + LATENCY + 2):
        await RisingEdge(dut.clk)
        dut.in_valid.value = int(i < len(xs))
        dut.x.value = mask(xs[i] if i < len(xs) else 0, W)

        await FallingEdge(dut.clk)
        if int(dut.out_valid.value):
            got.append((dut.y.value.to_signed(), int(dut.zero.value)))

    assert len(got) == len(golden), f"got {len(got)} results, expected {len(golden)}"
    for i, (g, ref) in enumerate(zip(got, golden)):
        assert g == ref, (f"sample {i}: x={xs[i]} = {xs[i] * 2.0 ** -F:.9f}  "
                          f"rtl=(y={g[0]}, zero={g[1]})  "
                          f"model=(y={ref[0]}, zero={ref[1]})  diff={g[0] - ref[0]}")
