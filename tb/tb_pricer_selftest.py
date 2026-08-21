"""cocotb testbench for pricer_selftest. The harness streams a ROM of vectors
through the pricer at one per clock and compares on chip, so the whole thing has
to be right before it is worth building: a self test that reports success
incorrectly is worse than no self test.

The negative control runs separately, because $readmemh loads the ROM once at
elaboration: corrupting the file mid simulation would change nothing."""

import os

import cocotb
from cocotb.clock import Clock
from cocotb.triggers import ClockCycles, RisingEdge

N = int(os.environ.get("N", "256"))
MEM_DIR = os.environ.get("MEM_DIR", ".")


async def run(dut, limit=20000):
    cocotb.start_soon(Clock(dut.clk, 10, "ns").start())
    dut.rst.value = 1
    dut.start.value = 0
    await ClockCycles(dut.clk, 5)
    dut.rst.value = 0
    await ClockCycles(dut.clk, 2)

    dut.start.value = 1
    await RisingEdge(dut.clk)
    dut.start.value = 0

    for _ in range(limit):
        await RisingEdge(dut.clk)
        if int(dut.done.value):
            break
    else:
        raise TimeoutError(
            f"done never asserted, seen={int(dut.seen.value)} of {N}")

    await ClockCycles(dut.clk, 2)
    return {
        "seen": int(dut.seen.value),
        "mismatches": int(dut.mismatches.value),
        "first_bad": int(dut.first_bad.value),
        "call": int(dut.first_bad_call.value),
        "put": int(dut.first_bad_put.value),
    }


@cocotb.test()
async def test_all_vectors_match(dut):
    """Every vector through the pipeline, compared on chip."""
    r = await run(dut)
    assert r["seen"] == N, f"only {r['seen']} of {N} results came out"
    assert r["mismatches"] == 0, (
        f"{r['mismatches']} mismatches, first at index {r['first_bad']}: "
        f"call={r['call']} put={r['put']}")
    dut._log.info(f"{r['seen']} vectors streamed, 0 mismatches")


@cocotb.test()
async def test_throughput_is_one_per_clock(dut):
    """The whole point of the pipeline: results arrive back to back."""
    cocotb.start_soon(Clock(dut.clk, 10, "ns").start())
    dut.rst.value = 1
    dut.start.value = 0
    await ClockCycles(dut.clk, 5)
    dut.rst.value = 0
    await ClockCycles(dut.clk, 2)
    dut.start.value = 1
    await RisingEdge(dut.clk)
    dut.start.value = 0

    first = last = None
    cycle = 0
    prev = 0
    for _ in range(20000):
        await RisingEdge(dut.clk)
        cycle += 1
        now = int(dut.seen.value)
        if now != prev:
            if first is None:
                first = cycle
            last = cycle
            prev = now
        if int(dut.done.value):
            break

    span = last - first + 1
    assert span == N, (
        f"{N} results spread over {span} clocks, expected {N} back to back")
    dut._log.info(f"{N} results in {span} clocks, one per clock")
