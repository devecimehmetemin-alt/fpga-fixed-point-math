"""cocotb testbench for exp_unit. drive a sweep, compare every output to the
golden model bit for bit."""

import os

import numpy as np

import cocotb
from cocotb.clock import Clock
from cocotb.triggers import ClockCycles, FallingEdge, RisingEdge

import exp_unit_model as exp
from fixedpoint import Fmt

X_W, X_F = int(os.environ["X_W"]), int(os.environ["X_F"])
L_W, L_F = int(os.environ["L_W"]), int(os.environ["L_F"])
T_W, T_F = int(os.environ["T_W"]), int(os.environ["T_F"])
C_W, C_F = int(os.environ["C_W"]), int(os.environ["C_F"])
W_W, W_F = int(os.environ["W_W"]), int(os.environ["W_F"])
Y_W, Y_F = int(os.environ["Y_W"]), int(os.environ["Y_F"])

MAX_SHIFT = 1
LATENCY = 12

def mask(value, width):
    # A Python int to wire
    return int(value) & ((1 << width) - 1)

@cocotb.test()
async def test_sweep(dut):
    FMT = exp.ExpFormats(X_W=X_W, X_F=X_F, L_W=L_W, L_F=L_F,
                                T_W=T_W, T_F=T_F, C_W=C_W, C_F=C_F,
                                W_W=W_W, W_F=W_F, Y_W=Y_W, Y_F=Y_F)

    X_MIN = -(2.0 ** (X_W - 1 - X_F))

    xs = np.concatenate([
        np.linspace(X_MIN, 0.0, 257),   
        exp.expo_boundaries(FMT),           
    ])
    stimulus = [int(c) for c in FMT.x.const(xs).code]
    golden   = [int(c) for c in exp.evaluate_fixed(xs, FMT)]

    cocotb.start_soon(Clock(dut.clk, 10 ,"ns").start())
    dut.rst.value = 1
    dut.in_valid.value = 0
    dut.x.value = 0

    await ClockCycles(dut.clk, 3)
    dut.rst.value = 0

    got = [] 
    for i in range(len(stimulus) + LATENCY + 2):
        await RisingEdge(dut.clk)
        dut.in_valid.value = int(i<len(stimulus))
        dut.x.value = mask(stimulus[i] if i < len(stimulus) else 0, X_W)

        await FallingEdge(dut.clk)
        if int(dut.out_valid.value):
            got.append(dut.y.value.to_unsigned())

    assert len(got) == len(golden), f"got {len(got)} results, expected {len(golden)}"
    for i, (g, e) in enumerate(zip(got, golden)):
        assert g == e, f"sample {i}: x={xs[i]:+.6f}  rtl={g}  model={e}  diff={g - e}"







