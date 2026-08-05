"""cocotb testbench for poly_eval. drive a sweep, compare every output to the
golden model bit for bit."""

import os

import numpy as np

import cocotb
from cocotb.clock import Clock
from cocotb.triggers import ClockCycles, FallingEdge, RisingEdge

import poly_eval_model as poly
from fixedpoint import Fmt

# Widths come from the Makefile, which passes the same numbers to Verilator as -G parameters. 
X_W, X_F = int(os.environ["X_W"]), int(os.environ["X_F"])
C_W, C_F = int(os.environ["C_W"]), int(os.environ["C_F"])
W_W, W_F = int(os.environ["W_W"]), int(os.environ["W_F"])

LATENCY = 7  # match the localparam in poly_eval.sv

X_FMT = Fmt(X_W, X_F, "half_up")
C_FMT = Fmt(C_W, C_F, "half_up")
W_FMT = Fmt(W_W, W_F, "half_up")

COEFFS = poly.exp2_coefficients()


def mask(value, width):
    # A Python int to wire value
    return int(value) & ((1 << width) - 1)


@cocotb.test()
async def test_sweep(dut):
    xs = np.linspace(-0.5, 0.5, 257)
    stimulus = [int(c) for c in X_FMT.const(xs).code]
    golden = [int(c) for c in poly.evaluate_estrin_fixed(COEFFS, xs, C_FMT, X_FMT, W_FMT).code]

    cocotb.start_soon(Clock(dut.clk, 10, "ns").start())
    dut.rst.value = 1
    dut.in_valid.value = 0
    dut.x.value = 0

    # dut.a is one 6*C_W vector with a0 in the low bits.
    word = 0
    for i, c in enumerate(poly.pad_coeffs(COEFFS)):
        word |= mask(C_FMT.const(c).code, C_W) << (i * C_W)
    dut.a.value = word
    await ClockCycles(dut.clk, 3)
    dut.rst.value = 0

    got = []
    for i in range(len(stimulus) + LATENCY + 2):
        # just after the rising edge
        await RisingEdge(dut.clk)
        dut.in_valid.value = int(i < len(stimulus))
        dut.x.value = mask(stimulus[i] if i < len(stimulus) else 0, X_W)

        # Sample at the falling edge
        await FallingEdge(dut.clk)
        if int(dut.out_valid.value):
            got.append(dut.p.value.to_signed())

    assert len(got) == len(golden), f"got {len(got)} results, expected {len(golden)}"
    for i, (g, e) in enumerate(zip(got, golden)):
        assert g == e, f"sample {i}: x={xs[i]:+.6f}  rtl={g}  model={e}  diff={g - e}"
