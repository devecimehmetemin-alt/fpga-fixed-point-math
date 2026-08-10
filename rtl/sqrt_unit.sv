// sqrt(a) by Newton-Raphson, for a = M * 2^e with M in [1,4) from lzc_norm.
// y_{n+1} = y_n(1.5 - 0.5*M*y_n*y_n)
// A midpoint seed table indexed by e[0] and the top IDX_W-1 bits of m starts
// it at SEED_BITS, so one iteration reaches 2*SEED_BITS.
// g1 = g0(1.5 - 0.5*M*y0*y0)
//    = g0(1 + 0.5*(1 - g0*y0))
//    = g0(1 + r)              r = 0.5*(1 - g0*y0), g0 = M*y0
//    = g0 + g0*r

module sqrt_unit #(
    parameter int A_W = 32, // input a width
    parameter int A_F = 28,
    parameter int Y_W = 24, // output y width
    parameter int Y_F = 21,
    parameter int IDX_W = 10, // seed table holds 2**IDX_W entries
    parameter int S_F = 17, // seed fractional bits
    parameter int M_W = 25, // bits of M the multiplier sees
    parameter int R_W = 18, // residual width
    parameter int E_MIN = 28, // declared exponent window
    parameter int E_MAX = 29,
    parameter int LATENCY = 12
)(
    input logic clk,
    input logic rst,
    input logic in_valid,
    input logic [A_W-1:0] a,
    output logic out_valid,
    output logic signed [Y_W-1:0] y
);

    localparam int E_W = $clog2(A_W) + 1; // matches lzc_norm
    localparam int LZC_LATENCY = 2;
    localparam int M_F = M_W - 3; // M is in [1,4): two integer bits, one sign

    localparam int NSEED = 1 << IDX_W;
    localparam int HSEED = 1 << (IDX_W - 1); // entries per parity half

    // g = M*y0 is sqrt(M), just over [1,2) because the seed can overshoot. The
    // full product feeds the C port; g_q feeds the multiplier ports and goes on
    // the wide one, because r measures how far g_q*y0 is from 1 while the
    // correction is applied to the full g -- a g_q error reaches the answer
    // undivided, unlike a seed error.
    localparam int PROD_W = M_W + S_F + 1; // M_q*y0
    localparam int G_F = M_F + S_F;
    localparam int GB_F = M_W - 3; // g reaches 2.001, so two integer bits
    localparam int SHIFT_G = G_F - GB_F;

    // as much headroom as the 48-bit accumulator allows on the C port
    localparam int SHIFT_C = 47 - PROD_W;

    // residual fractional bits
    localparam int R_F = G_F + SHIFT_C - GB_F;

    // product scale down to r's scale. The extra 1 is the 0.5 in the recurrence:
    // halving r once here costs nothing.
    localparam int SHIFT_RES = GB_F + S_F + 1 - R_F;

    localparam int G1_F = GB_F + R_F; // g0 + g0*r
    localparam int G1_W = PROD_W + SHIFT_C + 1;

    // Undoing the range reduction. The exponent is halved, so the window is too
    // and e[0] never reaches it: floor(e/2) already drops the odd bit.
    localparam int SHIFT_SPAN = (E_MAX >> 1) - (E_MIN >> 1);
    localparam int SHIFT_W = (SHIFT_SPAN < 1) ? 1 : $clog2(SHIFT_SPAN + 1);
    localparam int ACC_W = G1_W + SHIFT_SPAN;
    // The fixed shift is the one e=E_MIN needs, because a bigger e wants a
    // smaller shift here: sqrt grows with a where the reciprocal shrank.
    localparam int SHIFT_RIGHT = G1_F - Y_F - ((E_MIN >> 1) - A_F / 2);

    // 1.0 at the scale g_q*y0 comes out at
    localparam int PROD2_W = M_W + S_F + 2;
    localparam logic signed [PROD2_W-1:0] ONE = 1 <<< (GB_F + S_F);
    localparam logic signed [ACC_W-1:0] RND_Y = 1 <<< (SHIFT_RIGHT - 1);
    localparam logic signed [Y_W-1:0] Y_MAX = (1 <<< (Y_W - 1)) - 1;

    // Seed maker
    // Entry i holds the inverse square root of slice i's midpoint
    function automatic logic [NSEED-1:0][S_F-1:0] make_seed();
        logic [NSEED-1:0][S_F-1:0] t;
        real mid;
        int half;
        for (int i = 0; i < NSEED; i++) begin
            half = i >> (IDX_W - 1);
            mid = (1.0 + (real'(i & (HSEED - 1)) + 0.5) / real'(HSEED)) * (2.0 ** half);
            t[i] = S_F'($rtoi(1.0 / $sqrt(mid) * (2.0 ** S_F) + 0.5));
        end
        return t;
    endfunction

    localparam logic [NSEED-1:0][S_F-1:0] SEED = make_seed();

    // m's bits below m_q's slice never reach the multiplier.
    logic [A_W-1:0] m;
    logic signed [E_W-1:0] e;
    logic zero_lzc, lzc_valid;

    logic signed [M_W-1:0] m_q, M_q;
    logic [IDX_W-1:0] idx;
    logic signed [S_F:0] y0, y0_d [0:1];
    logic signed [PROD_W-1:0] prod, prod_d [0:4];
    logic signed [M_W-1:0] g_q, g_q_d [0:2];
    logic signed [PROD2_W-1:0] prod2, r_raw;
    logic signed [R_W-1:0] r;
    logic signed [M_W+R_W-1:0] prod3;
    logic signed [G1_W-1:0] g1;
    logic [SHIFT_W-1:0] shift_left;

    // exponent
    localparam int ED = LATENCY - LZC_LATENCY - 1;
    logic [SHIFT_W-1:0] shift_left_d [0:ED-1];
    logic [ED-1:0] sat_hi_d, sat_lo_d;
    logic sat_hi, sat_lo;
    logic [LATENCY-LZC_LATENCY-1:0] v;

    // clock 0-1
    lzc_norm #(
        .W (A_W),
        .LATENCY (LZC_LATENCY)
    ) u_lzc_norm (
        .clk (clk),
        .rst (rst),
        .in_valid (in_valid),
        .x (a),
        .out_valid (lzc_valid),
        .m (m),
        .e (e),
        .zero (zero_lzc)
    );

    always_ff @(posedge clk) begin
        // clock 2
        // idx samples e[0] and fractional bits of m
        // m_q samples the top M_W-2 bits of m, zero extended for the sign
        m_q <= signed'({2'b0, m[A_W-1 -: M_W-2]});
        idx <= {e[0], m[A_W-2 -: IDX_W-1]};

        // clock 3
        // an odd exponent moves M into [2,4), which is the spare integer bit.
        // The parity has to come from idx, not from e: e is lzc_norm's output
        // and by now holds the next sample's exponent.
        y0 <= signed'({1'b0, SEED[idx]});
        M_q <= m_q <<< idx[IDX_W-1];

        // clock 4
        // First multiply in Newton-Raphson, M*y0
        prod <= M_q * y0; // = sqrt(M)

        // clock 5
        // g down to the multiplier ports
        g_q <= M_W'(prod >>> SHIFT_G);

        // clock 6
        // g0*y0, which is M*y0*y0 and therefore near 1
        prod2 <= g_q * y0_d[1];

        // clock 7
        // 1 - g0*y0
        r_raw <= ONE - prod2;

        // clock 8
        // point alignment
        r <= R_W'(r_raw >>> SHIFT_RES);

        // clock 9-10
        // DSP3, g1 = g0 + g0*r, in form ax + b
        // then g0 + g0*r
        prod3 <= g_q_d[2] * r;
        g1 <= (G1_W'(prod_d[4]) <<< SHIFT_C) + G1_W'(prod3);

        // y0 is made at clock 3 and wanted at 6, g_q at 5 and wanted at 9,
        // prod at 4 and wanted at 10, therefore delay
        y0_d[0] <= y0;
        for (int i = 1; i < 2; i++) y0_d[i] <= y0_d[i-1];
        g_q_d[0] <= g_q;
        for (int i = 1; i < 3; i++) g_q_d[i] <= g_q_d[i-1];
        prod_d[0] <= prod;
        for (int i = 1; i < 5; i++) prod_d[i] <= prod_d[i-1];

        // clock 11
        // reconstruct
        // a = 0 gives e = 0, which is already below the window, so sqrt(0) = 0
        // falls out of the underflow branch and needs no flag of its own.
        if (sat_hi_d[ED-1]) begin
            y <= Y_MAX;
        end else if (sat_lo_d[ED-1]) begin
            y <= '0;
        end else begin
            y <= Y_W'((((ACC_W'(g1) <<< shift_left_d[ED-1]) + RND_Y) >>> SHIFT_RIGHT));
        end

        shift_left_d[0] <= shift_left;
        for (int i = 1; i < ED; i++) shift_left_d[i] <= shift_left_d[i-1];

        if (rst) begin
            v <= '0;
            sat_hi_d <= '0;
            sat_lo_d <= '0;
        end else begin
            v <= {v[LATENCY-LZC_LATENCY-2:0], lzc_valid};
            sat_hi_d <= {sat_hi_d[ED-2:0], sat_hi};
            sat_lo_d <= {sat_lo_d[ED-2:0], sat_lo};
        end
    end

    always_comb begin
        // sqrt is increasing, so the two ends swap against recip_unit: a small
        // a makes a small root, not a huge one.
        sat_hi = e > E_W'(E_MAX);
        sat_lo = e < E_W'(E_MIN);
        if (sat_lo) shift_left = '0;
        else if (sat_hi) shift_left = SHIFT_W'(SHIFT_SPAN);
        else shift_left = SHIFT_W'((int'(e) >>> 1) - (E_MIN >> 1));
    end

    assign out_valid = v[LATENCY-LZC_LATENCY-1];

endmodule
