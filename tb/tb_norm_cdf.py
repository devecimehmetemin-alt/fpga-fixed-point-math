"""cocotb testbench for norm_cdf. Drive the corner set plus a sweep, compare
every output to the golden model bit for bit, then check the two structural
properties that a bit-exact compare cannot see on its own."""

import os

import cocotb
from cocotb.clock import Clock
from cocotb.triggers import ClockCycles, FallingEdge, RisingEdge

import norm_cdf_model as nc

IN_W, IN_F = int(os.environ["IN_W"]), int(os.environ["IN_F"])
N_W = int(os.environ["N_W"])
LATENCY = int(os.environ["LATENCY"])

FMT = nc.NormCdfFormats(
    IN_W=IN_W, IN_F=IN_F, CLAMP=float(os.environ["CLAMP"]),
    P_W=int(os.environ["P_W"]), P_F=int(os.environ["P_F"]),
    T_W=int(os.environ["T_W"]), T_F=int(os.environ["T_F"]),
    C_W=int(os.environ["C_W"]), C_F=int(os.environ["C_F"]),
    W_W=int(os.environ["W_W"]), W_F=int(os.environ["W_F"]),
    U_W=int(os.environ["U_W"]), U_F=int(os.environ["U_F"]),
    PHI_W=int(os.environ["PHI_W"]), PHI_F=int(os.environ["PHI_F"]),
    N_W=N_W, N_F=int(os.environ["N_F"]),
    IDX_W=int(os.environ["IDX_W"]),
    EC_W=int(os.environ["EC_W"]), EC_F=int(os.environ["EC_F"]),
    EW_W=int(os.environ["EW_W"]), EW_F=int(os.environ["EW_F"]),
)

def mask(value, width):
    # A Python int to wire
    return int(value) & ((1 << width) - 1)


async def run(dut, stimulus):
    """Drive stimulus back to back, return one output per input."""
    cocotb.start_soon(Clock(dut.clk, 10, "ns").start())
    dut.rst.value = 1
    dut.in_valid.value = 0
    dut.x.value = 0

    await ClockCycles(dut.clk, 3)
    dut.rst.value = 0

    got = []
    for i in range(len(stimulus) + LATENCY + 2):
        await RisingEdge(dut.clk)
        dut.in_valid.value = int(i < len(stimulus))
        dut.x.value = mask(stimulus[i] if i < len(stimulus) else 0, IN_W)

        await FallingEdge(dut.clk)
        if int(dut.out_valid.value):
            got.append(dut.n.value.to_unsigned())
    return got


@cocotb.test()
async def test_corners_and_sweep(dut):
    # corners() straddles zero, the sign flip, both clamp edges, and
    # 1 + p|x| = 2 where recip_unit's exponent steps from E_MIN to E_MAX. That
    # last one is invisible to a uniform sweep. The sweep is layered on top so
    # the interior gets covered too.
    hi = FMT.CLAMP_CODE
    count = 1201
    sweep = [int(round(-hi + 2 * hi * i / (count - 1))) for i in range(count)]
    stimulus = sorted(set(nc.corners(FMT)) | set(sweep))

    golden = nc.evaluate_fixed(stimulus, FMT)
    got = await run(dut, stimulus)

    assert len(got) == len(golden), \
        f"got {len(got)} results, expected {len(golden)} -- LATENCY is wrong"

    for i, (g, e) in enumerate(zip(got, golden)):
        assert g == e, (
            f"sample {i}: x={stimulus[i] * 2.0 ** -IN_F:+.9f}"
            f"  rtl={g}  model={e}  diff={g - e}")

    # Structural, checked on the RTL's own outputs rather than inherited from
    # the model: a CDF that steps backwards is invisible to a bit-exact compare
    # if the model steps backwards too, and it is lethal to the Phi^-1
    # round-trip later.
    n_of = dict(zip(stimulus, got))
    drops = [(stimulus[i], got[i], got[i + 1])
             for i in range(len(got) - 1) if got[i + 1] < got[i]]
    assert not drops, f"{len(drops)} backward steps, first at {drops[0]}"

    # N(x) + N(-x) == 1 exactly, by construction. x = 0 is excluded: it is one
    # sample, not a pair, and it maps to the positive branch.
    pairs = [c for c in stimulus if c > 0 and -c in n_of]
    bad = [c for c in pairs if n_of[c] + n_of[-c] != FMT.ONE_N]
    assert not bad, (
        f"{len(bad)}/{len(pairs)} symmetry failures, first at x={bad[0]}: "
        f"{n_of[bad[0]]} + {n_of[-bad[0]]} != {FMT.ONE_N}")

    dut._log.info(
        f"{len(got)} samples bit-exact, {len(pairs)} symmetry pairs exact, "
        f"0 backward steps")
