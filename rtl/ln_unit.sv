// ln x = ln m + e*ln2, for x = m * 2^e with m in [1,2) from lzc_norm.
// m's MSB is the leading one, so z = m - 1 is a slice of m, which makes poly_eval
// fit ln(1+z) on [0,1) with no subtractor in front of it.


module ln_unit #(
    parameter int W = 32, // input x width
    parameter int F = 28,
    parameter int X_W = 18, // poly_eval argument width, sign bit + X_F
    parameter int X_F = 17,
    parameter int C_W = 18, // coefficient width
    parameter int C_F = 17,
    parameter int W_W = 18, // poly_eval working width
    parameter int W_F = 16,
    parameter int L_W = 18, // ln2 width
    parameter int Y_W = 24, // output y width, y has W_F fractional bits
    parameter int LATENCY = 11
)(
    input logic clk,
    input logic rst,
    input logic in_valid,
    input logic [W-1:0] x,
    output logic out_valid,
    output logic signed [Y_W-1:0] y,
    output logic zero
);

    localparam int E_W = $clog2(W) + 1;
    localparam int LZC_LATENCY = 2;
    localparam int POLY_LATENCY = 7;
    localparam int ACC_W = E_W + L_W + 1; // e*LN2, plus a carry bit
    localparam int ED = POLY_LATENCY; // e waits for clocks 2..8
    localparam int ZD = LATENCY - 3; // zero waits for clocks 2..9

    // error checks
    if (LATENCY != LZC_LATENCY + 1 + POLY_LATENCY + 1) begin : g_latency_check
        $error("LATENCY=%0d but the pipeline is %0d deep",
               LATENCY, LZC_LATENCY + 1 + POLY_LATENCY + 1);
    end
    if (X_W != X_F + 1) begin : g_arg_check
        $error("X_W=%0d must be X_F+1=%0d: z is in [0,1), one sign bit and no integer bits",
               X_W, X_F + 1);
    end
    if (X_F > W - 1) begin : g_slice_check
        $error("X_F=%0d exceeds the %0d fractional bits m carries", X_F, W - 1);
    end

    // LN2 carries W_F fractional bits
    localparam logic signed [L_W-1:0] LN2 = L_W'($rtoi(0.693147180560 * (2.0 ** W_F) + 0.5));

    // ln(1+z) on [0,1), minimax degree 5
    localparam logic signed [C_W-1:0] A0 = C_W'($rtoi(0.000009975033 * (2.0 ** C_F) + 0.5));
    localparam logic signed [C_W-1:0] A1 = C_W'($rtoi(0.999235483833 * (2.0 ** C_F) + 0.5));
    localparam logic signed [C_W-1:0] A2 = C_W'($rtoi(-0.490230723423 * (2.0 ** C_F) - 0.5));
    localparam logic signed [C_W-1:0] A3 = C_W'($rtoi(0.285272681091 * (2.0 ** C_F) + 0.5));
    localparam logic signed [C_W-1:0] A4 = C_W'($rtoi(-0.131581825089 * (2.0 ** C_F) - 0.5));
    localparam logic signed [C_W-1:0] A5 = C_W'($rtoi(0.030449004539 * (2.0 ** C_F) + 0.5));

    logic [W-1:0] m;
    logic signed [E_W-1:0] e;
    logic zero_lzc;

    logic signed [X_W-1:0] z;
    logic signed [5:0][C_W-1:0] a;
    logic signed [W_W-1:0] l;

    logic signed [E_W-1:0] e_d [0:ED-1];
    logic signed [ACC_W-1:0] prod;

    logic lzc_valid, z_valid, poly_valid;
    logic [ZD-1:0] zero_d;

    assign a = {A5, A4, A3, A2, A1, A0};

    // clock 0-1
    lzc_norm #(
        .W (W),
        .LATENCY (LZC_LATENCY)
    ) u_lzc_norm (
        .clk (clk),
        .rst (rst),
        .in_valid (in_valid),
        .x (x),
        .out_valid (lzc_valid),
        .m (m),
        .e (e),
        .zero (zero_lzc)
    );

    // clock 2-9
    always_ff @(posedge clk) begin
        // Dropping m's leading one is subtracting 1.
        z <= {1'b0, m[W-2 -: X_F]};

        // The exponent waits out the polynomial.
        e_d[0] <= e - E_W'(F);
        for (int i = 1; i < ED; i++) e_d[i] <= e_d[i-1];

        if (rst) begin
            z_valid <= 1'b0;
            zero_d <= '0;
        end else begin
            z_valid <= lzc_valid;
            zero_d <= {zero_d[ZD-2:0], zero_lzc};
        end
    end

    // clock 3-9
    poly_eval #(
        .X_W (X_W),
        .X_F (X_F),
        .C_W (C_W),
        .C_F (C_F),
        .W_W (W_W),
        .W_F (W_F)
    ) u_poly_eval (
        .clk (clk),
        .rst (rst),
        .in_valid (z_valid),
        .x (z),
        .a (a),
        .out_valid (poly_valid),
        .p (l)
    );

    always_ff @(posedge clk) begin
        // clock 9, e*ln(2)
        prod <= e_d[ED-1] * LN2;

        // clock 10, e*ln(2) + ln(m)
        y <= zero_d[ZD-1] ? '0 : Y_W'(prod + ACC_W'(l));
        zero <= zero_d[ZD-1];

        if (rst) out_valid <= 1'b0;
        else out_valid <= poly_valid;
    end

endmodule
