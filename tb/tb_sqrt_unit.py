"""cocotb testbench for sqrt_unit. drive a sweep, compare every output to the
golden model bit for bit."""

import os

import numpy as np

import cocotb
from cocotb.clock import Clock
from cocotb.triggers import ClockCycles, FallingEdge, RisingEdge

import sqrt_unit_model as sq

A_W, A_F = int(os.environ["A_W"]), int(os.environ["A_F"])
Y_W, Y_F = int(os.environ["Y_W"]), int(os.environ["Y_F"])
IDX_W, S_F = int(os.environ["IDX_W"]), int(os.environ["S_F"])
M_W, R_W = int(os.environ["M_W"]), int(os.environ["R_W"])
E_MIN, E_MAX = int(os.environ["E_MIN"]), int(os.environ["E_MAX"])

LATENCY = int(os.environ["LATENCY"])

FMT = sq.SqrtFormats(A_W=A_W, A_F=A_F, Y_W=Y_W, Y_F=Y_F, IDX_W=IDX_W,
                     S_F=S_F, M_W=M_W, R_W=R_W, E_MIN=E_MIN, E_MAX=E_MAX)


def mask(value, width):
    # A Python int to wire
    return int(value) & ((1 << width) - 1)


def stimulus():
    # corners lands on every seed slice boundary, both parities of e, and the
    # codes either side of the window
    xs = sq.corners(FMT)
    rng = np.random.default_rng(0)
    xs += [int(v) for v in rng.integers(FMT.lo_code, FMT.hi_code + 1, 512)]
    return xs


@cocotb.test()
async def test_sweep(dut):
    xs = stimulus()
    golden = sq.evaluate_fixed(xs, FMT)

    cocotb.start_soon(Clock(dut.clk, 10, "ns").start())
    dut.rst.value = 1
    dut.in_valid.value = 0
    dut.a.value = 0

    await ClockCycles(dut.clk, 3)
    dut.rst.value = 0

    got = []
    for i in range(len(xs) + LATENCY + 2):
        await RisingEdge(dut.clk)
        dut.in_valid.value = int(i < len(xs))
        dut.a.value = mask(xs[i] if i < len(xs) else 0, A_W)

        await FallingEdge(dut.clk)
        if int(dut.out_valid.value):
            got.append(dut.y.value.to_signed())

    assert len(got) == len(golden), f"got {len(got)} results, expected {len(golden)}"
    for i, (g, ref) in enumerate(zip(got, golden)):
        assert g == ref, (f"sample {i}: a={xs[i]} = {xs[i] * 2.0 ** -A_F:.9f}  "
                          f"rtl={g}  model={ref}  diff={g - ref}")
