// p(x) = a0 + a1*x + a2*x^2 + a3*x^3 + a4*x^4 + a5*x^5
// Factored as p = t0 + x2*(t1 + x2*t2), where
// x2 = x*x, t0 = a0 + a1*x, t1 = a2 + a3*x, t2 = a4 + a5*x

module poly_eval #(
    parameter int X_W = 18,
    parameter int X_F = 16,
    parameter int C_W = 18,
    parameter int C_F = 16,
    parameter int W_W = 18,
    parameter int W_F = 15
) (
    input logic clk,
    input logic rst,
    input logic in_valid,
    input logic signed [X_W-1:0] x,
    input logic signed [5:0][C_W-1:0] a,
    output logic out_valid,
    output logic signed [W_W-1:0] p
);

    localparam int LATENCY = 7;

    localparam int PX2_W = 2 * X_W; // x*x
    localparam int PAX_W = C_W + X_W; // a*x
    localparam int PW_W = 2 * W_W; // x2*t

    localparam int SUM1_W = PAX_W + 1; // a + a*x, one carry bit
    localparam int SUM2_W = PW_W + 1;
    localparam int ACC_W = (SUM1_W > SUM2_W) ? SUM1_W : SUM2_W;

    localparam int SHIFT_X2 = 2 * X_F - W_F;
    localparam int SHIFT_T = C_F + X_F - W_F;
    localparam int SHIFT_W = W_F;

    localparam logic signed [ACC_W-1:0] RND_X2 = ACC_W'(1) <<< (SHIFT_X2 - 1);
    localparam logic signed [ACC_W-1:0] RND_T = ACC_W'(1) <<< (SHIFT_T - 1);
    localparam logic signed [ACC_W-1:0] RND_W = ACC_W'(1) <<< (SHIFT_W - 1);

    logic signed [X_W-1:0] x_r;
    logic signed [5:0][C_W-1:0] a_r;

    (* use_dsp = "yes" *) logic signed [PX2_W-1:0] m_x2;
    (* use_dsp = "yes" *) logic signed [PAX_W-1:0] m1, m3, m5;
    logic signed [C_W-1:0] ad0, ad2, ad4;

    logic signed [W_W-1:0] x2, t0, t1, t2;

    (* use_dsp = "yes" *) logic signed [PW_W-1:0] m_w, m_p;
    logic signed [W_W-1:0] w;
    logic signed [W_W-1:0] t1_d;
    logic signed [W_W-1:0] x2_d1, x2_d2;
    logic signed [W_W-1:0] t0_d1, t0_d2, t0_d3;

    logic [LATENCY-1:0] v;

    always_ff @(posedge clk) begin
        // clock 0
        x_r <= x;
        a_r <= a;

        // clock 1
        m_x2 <= x_r * x_r;
        m1 <= $signed(a_r[1]) * x_r;
        m3 <= $signed(a_r[3]) * x_r;
        m5 <= $signed(a_r[5]) * x_r;
        ad0 <= a_r[0];
        ad2 <= a_r[2];
        ad4 <= a_r[4];

        // clock 2
        x2 <= W_W'((ACC_W'(m_x2) + RND_X2) >>> SHIFT_X2);
        t0 <= W_W'((((ACC_W'(ad0) <<< X_F) + ACC_W'(m1)) + RND_T) >>> SHIFT_T);
        t1 <= W_W'((((ACC_W'(ad2) <<< X_F) + ACC_W'(m3)) + RND_T) >>> SHIFT_T);
        t2 <= W_W'((((ACC_W'(ad4) <<< X_F) + ACC_W'(m5)) + RND_T) >>> SHIFT_T);

        // clock 3
        m_w <= x2 * t2;
        t1_d <= t1;
        x2_d1 <= x2;
        t0_d1 <= t0;

        // clock 4
        w <= W_W'((((ACC_W'(t1_d) <<< W_F) + ACC_W'(m_w)) + RND_W) >>> SHIFT_W);
        x2_d2 <= x2_d1;
        t0_d2 <= t0_d1;

        // clock 5
        m_p <= x2_d2 * w;
        t0_d3 <= t0_d2;

        // clock 6
        p <= W_W'((((ACC_W'(t0_d3) <<< W_F) + ACC_W'(m_p)) + RND_W) >>> SHIFT_W);

        if (rst) v <= '0;
        else v <= {v[LATENCY-2:0], in_valid};
    end

    assign out_valid = v[LATENCY-1];

endmodule
