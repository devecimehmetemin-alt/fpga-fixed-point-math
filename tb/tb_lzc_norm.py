"""cocotb testbench for lzc_norm. drive a sweep, compare every output to the
golden model bit for bit."""

import os

import numpy as np

import cocotb
from cocotb.clock import Clock
from cocotb.triggers import ClockCycles, FallingEdge, RisingEdge

import lzc_norm_model as lzc

W = int(os.environ["W"])
LATENCY = int(os.environ["LATENCY"])


def mask(value, width):
    # A Python int to wire
    return int(value) & ((1 << width) - 1)


@cocotb.test()
async def test_sweep(dut):
    xs = list(range(1 << W)) if W <= 12 else (
        lzc.corners(W) + [int(v) for v in np.random.default_rng(0).integers(0, 1 << W, 512)])
    golden = [(m, e, int(z)) for m, e, _, z in (lzc.normalize(x, W) for x in xs)]

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
            got.append((dut.m.value.to_unsigned(),
                        dut.e.value.to_signed(),
                        int(dut.zero.value)))

    assert len(got) == len(golden), f"got {len(got)} results, expected {len(golden)}"
    for i, (g, ref) in enumerate(zip(got, golden)):
        assert g == ref, f"sample {i}: x={xs[i]}  rtl={g}  model={ref}"
