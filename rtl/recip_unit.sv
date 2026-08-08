// 1/a by Newton-Raphson, for a = m * 2^e with m in [1,2) from lzc_norm.
// y_{n+1} = y_n(2 - m*y_n)
// A midpoint seed table indexed by the top IDX_W bits of m starts
// it at SEED_BITS, so one iteration reaches 2*SEED_BITS.
// y1 = y0(2 - m*y0)
//    = y0(1 + (1 - m*y0))
//    = y0(1 + r)              r = 1 - m*y0
//    = y0 + y0*r

module recip_unit #(
    parameter int A_W = 32, // input a width
    parameter int A_F = 28,
    parameter int Y_W = 24, // output y width
    parameter int Y_F = 22,
    parameter int IDX_W = 9, // seed table holds 2**IDX_W entries
    parameter int S_F = 17, // seed fractional bits
    parameter int M_W = 25, // bits of m the multiplier sees
    parameter int R_W = 18, // residual width
    parameter int E_MIN = 28, // declared exponent window
    parameter int E_MAX = 29,
    parameter int LATENCY = 10
)(
    input logic clk,
    input logic rst,
    input logic in_valid,
    input logic [A_W-1:0] a,
    output logic out_valid,
    output logic signed [Y_W-1:0] y,
    output logic ovf
);

    localparam int E_W = $clog2(A_W) + 1; // matches lzc_norm
    localparam int LZC_LATENCY = 2;
    localparam int M_F = M_W - 2; // m is in [1,2): one integer bit, one sign

    localparam int SEED_BITS = IDX_W + 1;
    localparam int SHIFT_W = $clog2(E_MAX - E_MIN + 1); // output shifter select width

    // residual fractional bits
    localparam int R_F = SEED_BITS + R_W - 2;

    // product scale down to r's scale
    localparam int SHIFT_RES = M_F + S_F - R_F;

    localparam int NSEED = 1 << IDX_W;
    localparam int PROD_W = M_W + S_F + 1; // m_q*y0
    localparam int Y1_W = S_F + R_F + 2; // y0 + y0*r
    localparam int ACC_W = Y1_W + E_MAX - E_MIN;

    // Undoing the range reduction
    localparam int SHIFT_FIX = S_F + R_F - A_F - Y_F + E_MIN;
    localparam int SHIFT_RIGHT = SHIFT_FIX + E_MAX - E_MIN;

    // 1.0 at the scale m_q*y0 comes out at
    localparam logic signed [PROD_W-1:0] ONE = 1 <<< (M_F + S_F);
    localparam logic signed [ACC_W-1:0] RND_Y = 1 <<< (SHIFT_RIGHT - 1);
    localparam logic signed [Y_W-1:0] Y_MAX = (1 <<< (Y_W - 1)) - 1;

    if (SHIFT_FIX <= 0) begin : g_shift_check
        $error("SHIFT_FIX=%0d: the output shift runs the wrong way", SHIFT_FIX);
    end

    // Seed maker
    // Entry i holds the reciprocal of slice i's midpoint
    function automatic logic [NSEED-1:0][S_F-1:0] make_seed();
        logic [NSEED-1:0][S_F-1:0] t;
        real mid;
        for (int i = 0; i < NSEED; i++) begin
            mid = 1.0 + (real'(i) + 0.5) / real'(NSEED);
            t[i] = S_F'($rtoi(1.0 / mid * (2.0 ** S_F) + 0.5));
        end
        return t;
    endfunction

    localparam logic [NSEED-1:0][S_F-1:0] SEED = make_seed();

    // m's bits below m_q's slice never reach the multiplier.
    logic [A_W-1:0] m;
    logic signed [E_W-1:0] e;
    logic zero_lzc, lzc_valid;

    logic signed [M_W-1:0] m_q, m_q_d;
    logic [IDX_W-1:0] idx;
    logic signed [S_F:0] y0, y0_d [0:3];
    logic signed [PROD_W-1:0] prod, r_raw;
    logic signed [R_W-1:0] r;
    logic signed [S_F+R_W:0] prod2;
    logic signed [Y1_W-1:0] y1;
    logic [SHIFT_W-1:0] shift_left;

    // exponent
    localparam int ED = LATENCY - LZC_LATENCY - 1;
    logic [SHIFT_W-1:0] shift_left_d [0:ED-1];
    logic [ED-1:0] sat_hi_d, sat_lo_d, zero_d;
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
        // idx samples fractional bits of m
        // m_q samples the top M_W-1 bits of m, zero extended for the sign
        m_q <= signed'({1'b0, m[A_W-1 -: M_W-1]});
        idx <= m[A_W-2 -: IDX_W];

        // clock 3
        // m_q delayed a cycle so it meets y0 at the multiplier.
        y0 <= signed'({1'b0, SEED[idx]});
        m_q_d <= m_q;

        // clock 4
        // First multiply in Newton-Raphson, m*y0
        prod <= m_q_d * y0;

        // clock 5
        // 1- m*y0
        r_raw <= ONE - prod;

        // clock 6
        // point alignment
        r <= R_W'(r_raw >>> SHIFT_RES);

        // clock 7-8
        // DSP2, y1 = y0 + y0*r, in form ax + b
        // then y0 + y0*r
        prod2 <= y0_d[2] * r;
        y1 <= (Y1_W'(y0_d[3]) <<< R_F) + Y1_W'(prod2);

        // y0 is made at clock 3 and wanted at 7 and 8 therefore delay
        y0_d[0] <= y0;
        for (int i = 1; i < 4; i++) y0_d[i] <= y0_d[i-1];

        // clock 9
        // reconstruct
        ovf <= zero_d[ED-1] || sat_hi_d[ED-1];
        if (zero_d[ED-1] || sat_hi_d[ED-1]) begin
            y <= Y_MAX;
        end else if (sat_lo_d[ED-1]) begin
            y <= '0;
        end else begin
            y <= Y_W'((((ACC_W'(y1) <<< shift_left_d[ED-1]) + RND_Y) >>> SHIFT_RIGHT));
        end

        shift_left_d[0] <= shift_left;
        for (int i = 1; i < ED; i++) shift_left_d[i] <= shift_left_d[i-1];

        if (rst) begin
            v <= '0;
            zero_d <= '0;
            sat_hi_d <= '0;
            sat_lo_d <= '0;
        end else begin
            v <= {v[LATENCY-LZC_LATENCY-2:0], lzc_valid};
            zero_d <= {zero_d[ED-2:0], zero_lzc};
            sat_hi_d <= {sat_hi_d[ED-2:0], sat_hi};
            sat_lo_d <= {sat_lo_d[ED-2:0], sat_lo};
        end
    end

    always_comb begin
        sat_hi = e < E_W'(E_MIN);
        sat_lo = e > E_W'(E_MAX);
        if (sat_hi) shift_left = SHIFT_W'(E_MAX - E_MIN);
        else if (sat_lo) shift_left = '0;
        else shift_left = SHIFT_W'(E_MAX - e);
    end

    assign out_valid = v[LATENCY-LZC_LATENCY-1];

endmodule
