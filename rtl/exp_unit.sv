module exp_unit #(
    parameter int X_W = 18, // input x width
    parameter int X_F = 17, // input x fractional width
    parameter int L_W = 18, // log2(e) width
    parameter int L_F = 16, 
    parameter int T_W = 20, // argument t width
    parameter int T_F = 18, 
    parameter int C_W = 18, // coefficient width
    parameter int C_F = 16, 
    parameter int W_W = 18, // significand width
    parameter int W_F = 16, 
    parameter int Y_W = 32, // output y width
    parameter int Y_F = 31, 
    parameter int MAX_SHIFT = 1, // largest right shift, set by x range
    parameter int LATENCY = 12 
)(
    input logic clk,
    input logic rst,
    input logic in_valid,
    input logic signed [X_W-1:0] x,
    output logic out_valid,
    output logic [Y_W-1:0] y
);

    // integer bits of t = width of expo
    localparam int EXPO_W = T_W - T_F;

    // surplus frac bits in x*LOG2E
    localparam int SH_T = X_F + L_F - T_F;

    // frac bits gained, signif -> y
    localparam int SH_Y = Y_F - W_F;

    // bits in shift, one mux level each
    localparam int SHIFT_STAGES = $clog2(MAX_SHIFT + 1);

    // width of x*LOG2E, +1 carry for +RND_BIAS
    localparam int ACC_W = X_W + L_W + 1;

    // round half up and +1/2
    localparam logic signed [ACC_W-1:0] RND_BIAS = (ACC_W'(1) <<< (SH_T-1)) + (ACC_W'(1) <<< (SH_T+T_F-1));

    // log2(e) at L_F frac bits
    localparam logic signed [L_W-1:0] LOG2E = L_W'($rtoi(1.4426950408889634 * (2.0 ** L_F) + 0.5));

    // must match poly_eval's LATENCY
    localparam int POLY_LATENCY = 7;

    // 3 reduction stages + the polynomial + 2 shifter stages
    if (LATENCY != 3 + POLY_LATENCY + 2) begin : g_latency_check
        $error("LATENCY=%0d but the pipeline is %0d deep", LATENCY, 3 + POLY_LATENCY + 2);
    end

    // minimax fit of 2^f on [-1/2, 1/2], converted from base 10 to fixed point
    localparam logic signed [C_W-1:0] A0 = C_W'($rtoi(1.000000075495 * (2.0 ** C_F) + 0.5));
    localparam logic signed [C_W-1:0] A1 = C_W'($rtoi(0.693147206711 * (2.0 ** C_F) + 0.5));
    localparam logic signed [C_W-1:0] A2 = C_W'($rtoi(0.240221073558 * (2.0 ** C_F) + 0.5));
    localparam logic signed [C_W-1:0] A3 = C_W'($rtoi(0.055503272142 * (2.0 ** C_F) + 0.5));
    localparam logic signed [C_W-1:0] A4 = C_W'($rtoi(0.009676037098 * (2.0 ** C_F) + 0.5));
    localparam logic signed [C_W-1:0] A5 = C_W'($rtoi(0.001340043216 * (2.0 ** C_F) + 0.5));

    logic signed [EXPO_W-1:0] expo;
    logic signed [T_F-1:0] frac;
    logic signed [X_W-1:0] x_r;
    logic signed [X_W+L_W-1:0] t_prod;
    logic signed [T_W-1:0] t_biased;
    logic signed [5:0][C_W-1:0] a;
    logic signed [W_W-1:0] signif;
    logic [W_W+SH_Y:0] signif_algn;
    logic [W_W+SH_Y:0] sh1;
    logic [SHIFT_STAGES-1:0] rshift;
    logic [POLY_LATENCY-1:0][SHIFT_STAGES-1:0] rshift_d;
    logic [2:0] v_red;
    logic poly_valid, v_shift;

    poly_eval #(
        .X_W(T_F),
        .X_F(T_F),
        .C_W(C_W),
        .C_F(C_F),
        .W_W(W_W),
        .W_F(W_F)
        ) u_poly_eval (
        .clk,
        .rst,
        .in_valid(v_red[2]),
        .x(frac),
        .a,
        .out_valid(poly_valid),
        .p(signif)
    );

    assign a = {A5, A4, A3, A2, A1, A0};

    assign expo = t_biased[T_W-1 : T_F];
    assign frac = { ~t_biased[T_F-1], t_biased[T_F-2:0] };

    // position signif to the output's binary pointer
    assign signif_algn = (W_W+SH_Y+1)'($unsigned(signif)) << (SH_Y+1);

    assign rshift = rshift_d[POLY_LATENCY-1];

    always_ff @(posedge clk) begin

        // clock 0
        x_r <= x;
        v_red[0] <= in_valid;

        // clock 1
        t_prod <= x_r * LOG2E;  //log2y = xlog2e
        v_red[1] <= v_red[0];

        // clock 2
        t_biased <= T_W'((t_prod + RND_BIAS) >>> SH_T); //
        v_red[2] <= v_red[1];

        // clock 3-9, expo waits alongside poly_eval
        rshift_d[0] <= SHIFT_STAGES'(-expo);
        for (int i = 1; i < POLY_LATENCY; i++) rshift_d[i] <= rshift_d[i-1];

        // clock 10
        v_shift <= poly_valid;
        sh1 <= signif_algn >> rshift;

        // clock 11
        out_valid <= v_shift;
        y <= Y_W'((sh1 + 1) >> 1);

        if (rst) begin
            v_red <= '0;
            v_shift <= 1'b0;
            out_valid <= 1'b0;
        end

    end

endmodule
