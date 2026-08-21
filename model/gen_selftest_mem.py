import itertools
import sys

import numpy as np

import bs_pricer_fixed as fx
import bs_pricer_model as bs

FMT = fx.DEFAULT
N = 256

SK_W, R_W, SG_W, T_W, OUT_W = FMT.SK_W, FMT.R_W, FMT.SG_W, FMT.T_W, FMT.OUT_W
IN_W = 2 * SK_W + R_W + SG_W + T_W
EXP_W = 2 * OUT_W


def mask(value, width):
    return int(value) & ((1 << width) - 1)


def vectors(n=N, seed=20260820):
    corners = list(itertools.product(*[bs.RANGES[a] for a in bs.ARGS]))
    rng = np.random.default_rng(seed)
    drawn = list(zip(*bs.sample(max(n - len(corners), 0), rng)))
    return [fx.encode(*row, fmt=FMT) for row in (corners + drawn)][:n]


def pack_in(row):
    s, k, r, sg, tau = row
    return (mask(s, SK_W) << (SK_W + R_W + SG_W + T_W)
            | mask(k, SK_W) << (R_W + SG_W + T_W)
            | mask(r, R_W) << (SG_W + T_W)
            | mask(sg, SG_W) << T_W
            | mask(tau, T_W))


def pack_exp(pair):
    call, put = pair
    return mask(call, OUT_W) << OUT_W | mask(put, OUT_W)


def main(outdir, n=N):
    rows = vectors(n)
    expect = fx.evaluate_fixed(rows, fmt=FMT)

    with open(f"{outdir}/vectors_in.mem", "w") as f:
        for row in rows:
            f.write(f"{pack_in(row):0{IN_W // 4}x}\n")
    with open(f"{outdir}/vectors_exp.mem", "w") as f:
        for pair in expect:
            f.write(f"{pack_exp(pair):0{EXP_W // 4}x}\n")

    print(f"{len(rows)} vectors -> {outdir}/vectors_in.mem ({IN_W}b), "
          f"{outdir}/vectors_exp.mem ({EXP_W}b)")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "../tb",
         int(sys.argv[2]) if len(sys.argv) > 2 else N)
