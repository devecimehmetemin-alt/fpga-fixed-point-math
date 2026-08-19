// N(x), the standard normal CDF using Hastings Approximation.
//
// t = 1 / (1 + p*|x|),   p = 0.2316419
// N(|x|) = 1 - phi(|x|) * (b1*t + b2*t^2 + b3*t^3 + b4*t^4 + b5*t^5)
// phi(x) = e^(-x*x/2) / sqrt(2*pi)
//
// 1/sqrt(2*pi) folds into the coefficients, so max |b'| is 0.7266 and the
// polynomial needs no integer bit.
//
// N(-|x|) = 1 - N(|x|) makes the negative branch a mux on the last stage
// instead of a second subtract, and it makes N(x) + N(-x) = 1 exact.
//
// |x| is clamped to CLAMP before anything else. Past there N has already
// saturated, so the clamp is exact behaviour and not an approximation.
//
// Two branches. They meet at clock 21.
// Branch A is 1 + p|x| -> recip -> poly, and 20 clocks.
// Branch B is e^(-x*x/2), and 20 clocks.

module norm_cdf #(
    parameter int IN_W = 26, // input x width
    parameter int IN_F = 22,
    parameter real CLAMP = 6.0, // |x| past this saturates, exactly
    parameter int P_W = 18, // constant p
    parameter int P_F = 17,
    parameter int T_W = 25, // t reaches 1.0, so T_F is one short of T_W-1
    parameter int T_F = 23,
    parameter int C_W = 18, // coefficients
    parameter int C_F = 17,
    parameter int W_W = 25, // polynomial working format
    parameter int W_F = 24,
    parameter int U_W = 28, // u = -x*x/2
    parameter int U_F = 22,
    parameter int PHI_W = 32, // e^u
    parameter int PHI_F = 31,
    parameter int N_W = 26, // output
    parameter int N_F = 25,
    parameter int IDX_W = 10, // recip_unit's seed table
    parameter int EC_W = 23, // exp_unit's own coefficient format
    parameter int EC_F = 21,
    parameter int EW_W = 24, // exp_unit's own significand format
    parameter int EW_F = 22,
    parameter int LATENCY = 23
)(
    input logic clk,
    input logic rst,
    input logic in_valid,
    input logic signed [IN_W-1:0] x,
    output logic out_valid,
    output logic [N_W-1:0] n
);

    localparam real INV_SQRT_2PI = 0.398942280401432678;
    localparam real LOG2E = 1.442695040888963407;
    localparam real P = 0.2316419;

    // remove sign bit
    localparam int AX_W = IN_W - 1;
    localparam int AX_F = IN_F;

    // recip_unit's input format
    localparam int A_W = 32;
    localparam int A_F = 28;

    localparam int RECIP_LATENCY = 10;
    localparam int POLY_LATENCY = 7;
    localparam int EXP_LATENCY = 12;

    // delay timings
    localparam int LATENCY_EXP = 4 + RECIP_LATENCY + POLY_LATENCY + 2;
    localparam int DELAY_B = LATENCY_EXP - 2 - 3 - EXP_LATENCY;

    if (LATENCY != LATENCY_EXP) begin : g_latency_check
        $error("LATENCY=%0d but the pipeline is %0d deep", LATENCY, LATENCY_EXP);
    end

    localparam int SHIFT_ARG = P_F + AX_F - A_F;
    localparam int SHIFT_U = 2 * AX_F - U_F + 1; // the halving rides the shift
    localparam int SHIFT_N = PHI_F + W_F - N_F;

    // max shift that could be needed
    // largest exponent magnitude we might work with is -18
    // e^-18 = 2^(-18*log2(e)) = 2^-25.97
    localparam int MAX_SHIFT = $rtoi($ceil(0.5 * CLAMP * CLAMP * LOG2E));

    localparam int PA_W = P_W + AX_W; // p*|x|, and the folded constant
    localparam int SQ_W = 2 * AX_W;
    localparam int PN_W = PHI_W + W_W;
    localparam int SIGN_D = LATENCY - 2; // sign made at 1, wanted at 22

    localparam logic [AX_W-1:0] CLAMP_CODE = AX_W'($rtoi(CLAMP * (2.0 ** IN_F) + 0.5));
    localparam logic [P_W-1:0] P_CODE = P_W'($rtoi(P * (2.0 ** P_F) + 0.5));

    // the 1 + and the round combined
    localparam logic [PA_W-1:0] C_ARG = (PA_W'(1) << (A_F + SHIFT_ARG)) + (PA_W'(1) << (SHIFT_ARG - 1));

    localparam logic [SQ_W-1:0] RND_U = SQ_W'(1) << (SHIFT_U - 1);
    localparam logic [PN_W-1:0] RND_N = PN_W'(1) << (SHIFT_N - 1);
    localparam logic [N_W-1:0] ONE_N = N_W'(1) << N_F;

    // b_i / sqrt(2*pi)
    function automatic logic signed [5:0][C_W-1:0] make_coeffs();
        logic signed [5:0][C_W-1:0] c;
        real b [0:5];
        b[0] = 0.0;
        b[1] = 0.319381530;
        b[2] = -0.356563782;
        b[3] = 1.781477937;
        b[4] = -1.821255978;
        b[5] = 1.330274429;
        for (int i = 0; i < 6; i++) begin
            c[i] = C_W'($rtoi($floor(b[i] * INV_SQRT_2PI * (2.0 ** C_F) + 0.5)));
        end
        return c;
    endfunction

    localparam logic signed [5:0][C_W-1:0] COEFF = make_coeffs();

    logic [IN_W-1:0] x_mag;
    logic [AX_W-1:0] ax_next, ax, ax_delay [0:DELAY_B-1];
    logic sign;
    logic [SIGN_D-1:0] sign_d;

    logic [AX_W-1:0] ax_dsp;
    (* use_dsp = "yes" *) logic [PA_W-1:0] p_times_ax, arg_biased;
    logic [PA_W-1:0] arg_full;
    logic [A_W-1:0] arg;
    logic signed [T_W-1:0] t;
    logic signed [W_W-1:0] poly_q;

    (* use_dsp = "yes" *) logic [SQ_W-1:0] ax_squared;
    logic [SQ_W-1:0] half_sq;
    logic signed [U_W-1:0] u;
    logic [PHI_W-1:0] phi;

    logic [PN_W-1:0] phi_times_q, phi_q_full;
    logic [N_W-1:0] phi_q;

    logic [3:0] vld_arg;
    logic [DELAY_B+1:0] vld_exp;
    logic [1:0] vld_out;
    logic t_valid, poly_q_valid, phi_valid, arg_ovf;

    // clock 4-14
    recip_unit #(
        .A_W (A_W),
        .A_F (A_F),
        .Y_W (T_W),
        .Y_F (T_F),
        .IDX_W (IDX_W),
        .LATENCY (RECIP_LATENCY)
    ) u_recip_unit (
        .clk (clk),
        .rst (rst),
        .in_valid (vld_arg[3]),
        .a (arg),
        .out_valid (t_valid),
        .y (t),
        .ovf (arg_ovf) // arg >= 1 always, so neither zero nor saturation
    );

    // clock 14-21
    poly_eval #(
        .X_W (T_W),
        .X_F (T_F),
        .C_W (C_W),
        .C_F (C_F),
        .W_W (W_W),
        .W_F (W_F)
    ) u_poly_eval (
        .clk (clk),
        .rst (rst),
        .in_valid (t_valid),
        .x (t),
        .a (COEFF),
        .out_valid (poly_q_valid),
        .p (poly_q)
    );

    // clock 9-21
    exp_unit #(
        .X_W (U_W),
        .X_F (U_F),
        .L_W (24),
        .L_F (22),
        .T_W (U_F + 6),
        .T_F (U_F),
        .C_W (EC_W),
        .C_F (EC_F),
        .W_W (EW_W),
        .W_F (EW_F),
        .Y_W (PHI_W),
        .Y_F (PHI_F),
        .MAX_SHIFT (MAX_SHIFT),
        .LATENCY (EXP_LATENCY)
    ) u_exp_unit (
        .clk (clk),
        .rst (rst),
        .in_valid (vld_exp[DELAY_B+1]),
        .x (u),
        .out_valid (phi_valid),
        .y (phi)
    );

    // arg = 1 + p*|x|, rounded and added 1 inside arg_biased
    assign arg_full = arg_biased >> SHIFT_ARG;
    assign arg = A_W'(arg_full);

    always_comb begin
        // magnitude at the input width
        x_mag = x[IN_W-1] ? $unsigned(-x) : $unsigned(x);

        // clamp at needed width
        ax_next = (x_mag > IN_W'(CLAMP_CODE)) ? CLAMP_CODE : AX_W'(x_mag);

        // round, halve and make at the product's width (SHIFT_U has plus 1 for divide by 2)
        half_sq = (ax_squared + RND_U) >> SHIFT_U;

        // round then cast
        phi_q_full = (phi_times_q + RND_N) >> SHIFT_N;
        phi_q = N_W'(phi_q_full);
    end

    always_ff @(posedge clk) begin
        // clock 1
        sign <= x[IN_W-1];
        ax <= ax_next;

        // clock 2-4
        // the DSP's A, M and P registers
        ax_dsp <= ax;
        p_times_ax <= P_CODE * ax_dsp;
        arg_biased <= p_times_ax + C_ARG;

        // clock 2-7
        // branch B waiting
        ax_delay[0] <= ax;
        for (int i = 1; i < DELAY_B; i++) ax_delay[i] <= ax_delay[i-1];

        // clock 8
        ax_squared <= ax_delay[DELAY_B-1] * ax_delay[DELAY_B-1];

        // clock 9
        // -x^2 / 2, negate last so the shift is unsigned
        u <= -signed'(U_W'(half_sq));

        // clock 22
        // phi * polynomial output
        phi_times_q <= phi * $unsigned(poly_q);

        // clock 23
        // 1 - phi*q' for x >= 0, phi*q' for x < 0
        n <= sign_d[SIGN_D-1] ? phi_q : ONE_N - phi_q;

        sign_d <= {sign_d[SIGN_D-2:0], sign};

        if (rst) begin
            vld_arg <= '0;
            vld_exp <= '0;
            vld_out <= '0;
        end else begin
            vld_arg <= {vld_arg[2:0], in_valid};
            vld_exp <= {vld_exp[DELAY_B:0], vld_arg[0]};
            vld_out <= {vld_out[0], poly_q_valid};
        end
    end

    assign out_valid = vld_out[1];

endmodule
