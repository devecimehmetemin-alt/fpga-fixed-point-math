// Black-Scholes European call and put, closed form
//
// v  = sigma*sqrt(T)
// d1 = [ln(S/K) + (r + sigma^2/2)*T] / v
// d2 = d1 - v
// C  = S*N(d1) - K*e^(-rT)*N(d2)
// P  = K*e^(-rT)*N(-d2) - S*N(-d1)
//
// Every elementary function is a unit from this library: recip_unit twice,
// sqrt_unit, ln_unit, exp_unit, norm_cdf twice. What is left is scheduling and
// format alignment.
//
// No degenerate bypass. v >= 3.16e-04 over the box, and recip_unit saturating
// below its exponent window drives |d| up rather than to a nan, so the clamp
// delivers discounted intrinsic.
//
// Nine delay lines, none reset, so they infer as SRLs.
//
// Four branches, meeting in two joins. One result per clock, 53 of fill.
//   K:    1/K at 10, x at 11, ln(x) at 22
//   T:    sqrt(T) at 12, v at 13, 1/v at 23
//   R:    sigma^2 at 1, r + sigma^2/2 at 2, *T at 3
//   disc: r*T at 1, -rT at 2, e^(-rT) at 14
//   join 1 at 23: num = ln(x) + hvT, meeting 1/v
//   join 2 at 24: num*(1/v), d1 and d2 at 25, both clamps at 26
//   N(d1), N(d2) at 49, combine at 50-51, scale by K at 52, floor at 53

module bs_pricer_top #(
    parameter int SK_W = 32, // S and K share a format
    parameter int SK_F = 22,
    parameter int T_W = 32, // maturity
    parameter int T_F = 28,
    parameter int SG_W = 24, // volatility
    parameter int SG_F = 22,
    parameter int R_W = 24, // rate, signed: negative rates are in the box
    parameter int R_F = 22,
    parameter int IK_W = 32, // 1/K
    parameter int IK_F = 32,
    parameter int X_W = 36, // x = S/K
    parameter int X_F = 30,
    parameter int LX_W = 28, // ln(x)
    parameter int LX_F = 24,
    parameter int ST_W = 28, // sqrt(T)
    parameter int ST_F = 25,
    parameter int V_W = 28, // v = sigma*sqrt(T)
    parameter int V_F = 25,
    parameter int IV_W = 32, // 1/v
    parameter int IV_F = 19,
    parameter int HV_W = 26, // r + sigma^2/2
    parameter int HV_F = 24,
    parameter int HT_W = 28, // (r + sigma^2/2)*T
    parameter int HT_F = 24,
    parameter int NM_W = 29, // ln(x) + hvT
    parameter int NM_F = 24,
    parameter int D_W = 40, // d1, d2 before the clamp
    parameter int IN_W = 26, // norm_cdf's input, d after the clamp
    parameter int IN_F = 22,
    parameter real CLAMP = 6.0,
    parameter int N_W = 26, // norm_cdf's output
    parameter int N_F = 25,
    parameter int U_W = 28, // u = -rT
    parameter int U_F = 22,
    parameter int DC_W = 32, // disc = e^u
    parameter int DC_F = 31,
    parameter int CN_W = 32, // C/K
    parameter int CN_F = 26,
    parameter int PN_W = 28, // P/K
    parameter int PN_F = 26,
    parameter int OUT_W = 40, // C and P
    parameter int OUT_F = 22,
    parameter int K_IDX = 12, // seed table for 1/K
    parameter int V_IDX = 10, // seed table for 1/v
    parameter int EC_W = 23, // exp_unit's coefficient format, shared
    parameter int EC_F = 21,
    parameter int EW_W = 24, // exp_unit's significand format, shared
    parameter int EW_F = 22,
    parameter int LATENCY = 53
)(
    input logic clk,
    input logic rst,
    input logic in_valid,
    input logic [SK_W-1:0] s,
    input logic [SK_W-1:0] k,
    input logic signed [R_W-1:0] r,
    input logic [SG_W-1:0] sigma,
    input logic [T_W-1:0] tau,
    output logic out_valid,
    output logic signed [OUT_W-1:0] call,
    output logic signed [OUT_W-1:0] put
);

    localparam int RECIP_LATENCY = 10;
    localparam int SQRT_LATENCY = 12;
    localparam int LN_LATENCY = 11;
    localparam int EXP_LATENCY = 12;
    localparam int NCDF_LATENCY = 23;

    // clocks at which a signal lands
    localparam int C_IK = RECIP_LATENCY;
    localparam int C_X = C_IK + 1;
    localparam int C_LNX = C_X + LN_LATENCY;

    localparam int C_ST = SQRT_LATENCY;
    localparam int C_V = C_ST + 1;
    localparam int C_IV = C_V + RECIP_LATENCY;

    // Each step is registered, so the next calculation happens on the following clock.
    localparam int C_SG2 = 1;
    localparam int C_HV = C_SG2 + 1;
    localparam int C_HT = C_HV + 1;

    // Build -rT in two clocks, then pass it through the exponential unit.
    // Over the supported input range, exp_unit only needs a maximum shift of 1.
    localparam int C_U = 2;
    localparam int C_DC = C_U + EXP_LATENCY;
    localparam int MAX_SHIFT = 1;

    localparam int C_NUM = ((C_LNX > C_HT) ? C_LNX : C_HT) + 1;
    localparam int C_PD = ((C_NUM > C_IV) ? C_NUM : C_IV) + 1;

    // Spread the d1/d2 work across separate clocks to keep the critical path short.
    localparam int C_DW = C_PD + 1;
    localparam int C_D = C_DW + 1;
    localparam int C_ND = C_D + NCDF_LATENCY;
    localparam int C_PROD = C_ND + 1;
    localparam int C_CN = C_PROD + 1;
    localparam int C_PK = C_CN + 1;
    localparam int C_OUT = C_PK + 1;

    // Delay each signal until the clock where it is needed.
    localparam int DELAY_S = C_X - 1;
    localparam int DELAY_SG = C_V - 1;
    localparam int DELAY_TAU = C_HT - 1;
    localparam int DELAY_R = C_HV - 1;
    localparam int DELAY_HT = C_NUM - 1 - C_HT;

    // Convert v to the d1/d2 format early, then delay it until d2 is calculated.
    localparam int C_V2D = C_V + 1;
    localparam int DELAY_V2D = C_DW - 1 - C_V2D;

    localparam int DELAY_X = C_PROD - 1 - C_X;
    localparam int DELAY_DC = C_PROD - 1 - C_DC;
    localparam int DELAY_K = C_PK - 1;

    if (LATENCY != C_OUT) begin : g_latency_check
        $error("LATENCY=%0d but the pipeline is %0d deep", LATENCY, C_OUT);
    end
    if (C_NUM != C_IV) begin : g_join_check
        // num and 1/v feed the same multiply, should land at same clock
        $error("num lands at %0d and 1/v at %0d", C_NUM, C_IV);
    end
    if (HV_F < R_F) begin : g_rate_scale_check
        $error("HV_F=%0d is below R_F=%0d, so r cannot align by shifting left",
               HV_F, R_F);
    end

    // shift calculation for multiply ops
    localparam int SHIFT_X = SK_F + IK_F - X_F;
    localparam int SHIFT_V = SG_F + ST_F - V_F;
    localparam int SHIFT_HV = 2 * SG_F + 1 - HV_F;
    localparam int SHIFT_R = HV_F - R_F;
    localparam int SHIFT_HT = HV_F + T_F - HT_F;
    localparam int SHIFT_NM = LX_F - NM_F;
    localparam int SHIFT_D = NM_F + IV_F - IN_F;
    localparam int SHIFT_V2D = V_F - IN_F;
    localparam int V2D_W = V_W - SHIFT_V2D;
    localparam int SHIFT_U = R_F + T_F - U_F;
    localparam int SHIFT_XN = X_F + N_F - CN_F;
    localparam int SHIFT_DN = DC_F + N_F - CN_F;
    localparam int SHIFT_PN = DC_F + N_F - PN_F;
    localparam int SHIFT_XP = X_F + N_F - PN_F;
    localparam int SHIFT_C = SK_F + CN_F - OUT_F;
    localparam int SHIFT_P = SK_F + PN_F - OUT_F;

    // Use a minimum shift of 1 to avoid a negative shift
    // when calculating the rounding constant.
    localparam int SHIFT_NM_S = (SHIFT_NM > 0) ? SHIFT_NM : 1;

    // full bit widths of multiplations
    localparam int PX_W = SK_W + IK_W;
    localparam int PV_W = SG_W + ST_W;
    localparam int PSG_W = 2 * SG_W;
    localparam int PHT_W = HV_W + T_W + 1;
    localparam int PD_W = NM_W + IV_W;
    localparam int PU_W = R_W + T_W + 1;
    localparam int PXN_W = X_W + N_W;
    localparam int PDN_W = DC_W + N_W;
    localparam int PC_W = SK_W + CN_W + 1;
    localparam int PP_W = SK_W + PN_W + 1;

    // Add one extra bit so the call/put subtractions have room for a signed result.
    localparam int CS_W = CN_W + 1;
    localparam int PS_W = PN_W + 1;

    // add half for rounding
    localparam logic [PX_W-1:0] RND_X = PX_W'(1) << (SHIFT_X - 1);
    localparam logic [PV_W-1:0] RND_V = PV_W'(1) << (SHIFT_V - 1);
    localparam logic [PSG_W-1:0] RND_HV = PSG_W'(1) << (SHIFT_HV - 1);
    localparam logic signed [PHT_W-1:0] RND_HT = PHT_W'(1) << (SHIFT_HT - 1);
    localparam logic signed [LX_W-1:0] RND_NM = (SHIFT_NM > 0) ? (LX_W'(1) << (SHIFT_NM_S - 1)) : '0;
    localparam logic signed [PD_W-1:0] RND_D = PD_W'(1) << (SHIFT_D - 1);
    localparam logic [V_W-1:0] RND_V2D = V_W'(1) << (SHIFT_V2D - 1);
    localparam logic signed [PU_W-1:0] RND_U = PU_W'(1) << (SHIFT_U - 1);
    localparam logic [PXN_W-1:0] RND_XN = PXN_W'(1) << (SHIFT_XN - 1);
    localparam logic [PDN_W-1:0] RND_DN = PDN_W'(1) << (SHIFT_DN - 1);
    localparam logic [PDN_W-1:0] RND_PN = PDN_W'(1) << (SHIFT_PN - 1);
    localparam logic [PXN_W-1:0] RND_XP = PXN_W'(1) << (SHIFT_XP - 1);
    localparam logic signed [PC_W-1:0] RND_C = PC_W'(1) << (SHIFT_C - 1);
    localparam logic signed [PP_W-1:0] RND_P = PP_W'(1) << (SHIFT_P - 1);

    localparam logic signed [D_W-1:0] CLAMP_CODE = D_W'($rtoi(CLAMP * (2.0 ** IN_F) + 0.5));
    localparam logic [N_W-1:0] ONE_N = N_W'(1) << N_F;

    logic signed [IK_W-1:0] inv_k;
    logic [PX_W-1:0] p_x;
    logic [X_W-1:0] x;
    logic signed [LX_W-1:0] lnx;

    logic signed [ST_W-1:0] sq_t;
    logic [PV_W-1:0] p_v;
    logic [V_W-1:0] v;
    logic signed [IV_W-1:0] inv_v;

    logic [PSG_W-1:0] p_sg;
    logic signed [HV_W-1:0] hv_sq, hv;
    logic signed [PHT_W-1:0] p_ht;
    logic signed [HT_W-1:0] hv_t;

    // temporary full width results after rounding but before narrowing
    logic [PX_W-1:0] x_full;
    logic [PV_W-1:0] v_full;
    logic [PSG_W-1:0] sq_full;
    logic signed [PHT_W-1:0] ht_full;
    logic signed [LX_W-1:0] nm_full;
    logic signed [PD_W-1:0] d_full;
    logic signed [PU_W-1:0] u_full;

    // fractionally formatted
    logic [V2D_W-1:0] v2d;
    logic [N_W-1:0] cd1, cd2;
    logic signed [NM_W-1:0] num;
    logic signed [PD_W-1:0] p_d;
    logic signed [D_W-1:0] d1_wide, d2_wide, d1_sat, d2_sat;
    logic signed [IN_W-1:0] d1, d2;
    logic [N_W-1:0] nd1, nd2;

    logic signed [PU_W-1:0] p_u;
    logic signed [U_W-1:0] u;
    logic [DC_W-1:0] disc;

    logic [PXN_W-1:0] p_xn, p_xp;
    logic [PDN_W-1:0] p_dn, p_pn;
    logic [PXN_W-1:0] xn_full, xp_full;
    logic [PDN_W-1:0] dn_full, pn_full;
    logic signed [CS_W-1:0] call_s;
    logic signed [PS_W-1:0] put_s;
    logic signed [CN_W-1:0] call_n;
    logic signed [PN_W-1:0] put_n;

    logic signed [PC_W-1:0] p_c, call_raw;
    logic signed [PP_W-1:0] p_p, put_raw;

    // no reset, so they infer as SRLs
    logic [SK_W-1:0] s_d [0:DELAY_S-1];
    logic [SG_W-1:0] sg_d [0:DELAY_SG-1];
    logic [T_W-1:0] tau_d [0:DELAY_TAU-1];
    logic signed [R_W-1:0] r_d [0:DELAY_R-1];
    logic signed [HT_W-1:0] ht_d [0:DELAY_HT-1];
    logic [V2D_W-1:0] v2d_d [0:DELAY_V2D-1];
    logic [X_W-1:0] x_d [0:DELAY_X-1];
    logic [DC_W-1:0] dc_d [0:DELAY_DC-1];
    logic [SK_W-1:0] k_d [0:DELAY_K-1];

    logic [1:0] vld_u;
    logic vld_x, vld_v, vld_num, vld_pd, vld_dw, vld_d;
    logic ik_valid, st_valid, iv_valid, lnx_valid, dc_valid;
    logic nd1_valid, nd2_valid;
    logic [3:0] vld_back;
    logic ik_ovf, iv_ovf, lnx_zero;

    // clock 0-10
    recip_unit #(
        .A_W (SK_W),
        .A_F (SK_F),
        .Y_W (IK_W),
        .Y_F (IK_F),
        .IDX_W (K_IDX),
        .E_MIN (25),
        .E_MAX (30),
        .LATENCY (RECIP_LATENCY)
    ) u_recip_k (
        .clk (clk),
        .rst (rst),
        .in_valid (in_valid),
        .a (k),
        .out_valid (ik_valid),
        .y (inv_k),
        .ovf (ik_ovf) // K >= 10, so 1/K never saturates
    );

    // clock 11-22
    ln_unit #(
        .W (X_W),
        .F (X_F),
        .X_W (25),
        .X_F (24),
        .C_W (25),
        .C_F (24),
        .W_W (25),
        .W_F (LX_F),
        .L_W (26),
        .Y_W (LX_W),
        .LATENCY (LN_LATENCY)
    ) u_ln_unit (
        .clk (clk),
        .rst (rst),
        .in_valid (vld_x),
        .x (x),
        .out_valid (lnx_valid),
        .y (lnx),
        .zero (lnx_zero)
    );

    // clock 0-12
    sqrt_unit #(
        .A_W (T_W),
        .A_F (T_F),
        .Y_W (ST_W),
        .Y_F (ST_F),
        .IDX_W (10),
        .E_MIN (18),
        .E_MAX (30),
        .LATENCY (SQRT_LATENCY)
    ) u_sqrt_t (
        .clk (clk),
        .rst (rst),
        .in_valid (in_valid),
        .a (tau),
        .out_valid (st_valid),
        .y (sq_t)
    );

    // clock 13-23
    recip_unit #(
        .A_W (IV_W),
        .A_F (V_F),
        .Y_W (IV_W),
        .Y_F (IV_F),
        .IDX_W (V_IDX),
        .E_MIN (13),
        .E_MAX (26),
        .LATENCY (RECIP_LATENCY)
    ) u_recip_v (
        .clk (clk),
        .rst (rst),
        .in_valid (vld_v),
        .a (IV_W'(v)),
        .out_valid (iv_valid),
        .y (inv_v),
        .ovf (iv_ovf) // v stays inside [E_MIN,E_MAX]
    );

    // clock 2-14
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
        .Y_W (DC_W),
        .Y_F (DC_F),
        .MAX_SHIFT (MAX_SHIFT),
        .LATENCY (EXP_LATENCY)
    ) u_exp_disc (
        .clk (clk),
        .rst (rst),
        .in_valid (vld_u[1]),
        .x (u),
        .out_valid (dc_valid),
        .y (disc)
    );

    // clock 26-49. Two instances
    norm_cdf #(
        .IN_W (IN_W),
        .IN_F (IN_F),
        .CLAMP (CLAMP),
        .N_W (N_W),
        .N_F (N_F),
        .EC_W (EC_W),
        .EC_F (EC_F),
        .EW_W (EW_W),
        .EW_F (EW_F),
        .LATENCY (NCDF_LATENCY)
    ) u_norm_cdf_1 (
        .clk (clk),
        .rst (rst),
        .in_valid (vld_d),
        .x (d1),
        .out_valid (nd1_valid),
        .n (nd1)
    );

    norm_cdf #(
        .IN_W (IN_W),
        .IN_F (IN_F),
        .CLAMP (CLAMP),
        .N_W (N_W),
        .N_F (N_F),
        .EC_W (EC_W),
        .EC_F (EC_F),
        .EW_W (EW_W),
        .EW_F (EW_F),
        .LATENCY (NCDF_LATENCY)
    ) u_norm_cdf_2 (
        .clk (clk),
        .rst (rst),
        .in_valid (vld_d),
        .x (d2),
        .out_valid (nd2_valid),
        .n (nd2)
    );

    always_comb begin
        // All round-half-up at the product's width. A size cast sets the width
        // its operand is evaluated at, so W'(p >> S) chops p before the shift.
        x_full = (p_x + RND_X) >> SHIFT_X;
        x = X_W'(x_full);

        v_full = (p_v + RND_V) >> SHIFT_V;
        v = V_W'(v_full);

        sq_full = (p_sg + RND_HV) >> SHIFT_HV;
        hv_sq = signed'(HV_W'(sq_full));

        ht_full = (p_ht + RND_HT) >>> SHIFT_HT;
        hv_t = HT_W'(ht_full);

        nm_full = (lnx + RND_NM) >>> SHIFT_NM;

        u_full = (p_u + RND_U) >>> SHIFT_U;

        d_full = (p_d + RND_D) >>> SHIFT_D;

        // Saturate, then narrow. d1 runs to +10756 in the money and d2 to -10756
        // out of it, so both ends need it.
        d1_sat = (d1_wide > CLAMP_CODE) ? CLAMP_CODE
               : (d1_wide < -CLAMP_CODE) ? -CLAMP_CODE : d1_wide;
        d2_sat = (d2_wide > CLAMP_CODE) ? CLAMP_CODE
               : (d2_wide < -CLAMP_CODE) ? -CLAMP_CODE : d2_wide;

        // 1 - N is exact: norm_cdf's output format has 1.0 in it by
        // construction, so the complement is a subtract and not a rounding.
        cd1 = ONE_N - nd1;
        cd2 = ONE_N - nd2;

        xn_full = (p_xn + RND_XN) >> SHIFT_XN;
        dn_full = (p_dn + RND_DN) >> SHIFT_DN;
        pn_full = (p_pn + RND_PN) >> SHIFT_PN;
        xp_full = (p_xp + RND_XP) >> SHIFT_XP;

        // Both terms are non-negative but the difference is not, so widen before
        // the subtract.
        call_s = signed'(CS_W'(xn_full)) - signed'(CS_W'(dn_full));
        put_s = signed'(PS_W'(pn_full)) - signed'(PS_W'(xp_full));

        call_raw = (p_c + RND_C) >>> SHIFT_C;
        put_raw = (p_p + RND_P) >>> SHIFT_P;
    end

    always_ff @(posedge clk) begin
        // clock 1
        // sigma ^ 2 and rT
        p_sg <= sigma * sigma;
        p_u <= signed'(PU_W'(r)) * signed'({1'b0, tau});

        // clock 2
        // r + sigma ^ 2 / 2, hv_sq comes from combinational logic
        hv <= (HV_W'(r_d[DELAY_R-1]) <<< SHIFT_R) + hv_sq;

        // clock 2
        // negate rT
        u <= -U_W'(u_full);

        // clock 3
        // (r + sigma ^ 2/2)*T
        p_ht <= hv * signed'({1'b0, tau_d[DELAY_TAU-1]});

        // clock 11
        // recip unit out 1/k is multiplied by S
        // rounded and shifted combinationally before entering ln unit
        p_x <= s_d[DELAY_S-1] * $unsigned(inv_k);

        // clock 13
        // sqrt(T) is calculated at clock 12
        // calculate sigma * sqrt(T), rounded and scaled combinationally
        p_v <= sg_d[DELAY_SG-1] * $unsigned(sq_t);

        // clock 14
        // convert to format used by d2
        v2d <= V2D_W'((v + RND_V2D) >> SHIFT_V2D);

        // clock 23
        // 1/(sigma * sqrt(T)) finished calculating
        // num is set to ln(S/K) + (r + sigma^2/2)*T
        num <= NM_W'(nm_full) + NM_W'(ht_d[DELAY_HT-1]);

        // clock 24
        // num / (sigma * sqrt(T))
        // rounded and scaled down by combinational logic
        p_d <= num * inv_v;

        // clock 25
        // register rescaled result
        // d2 = d1 - sigma * sqrt(T)
        d1_wide <= D_W'(d_full);
        d2_wide <= D_W'(d_full) - signed'(D_W'(v2d_d[DELAY_V2D-1]));

        // clock 26
        // d1 d2 clamped combinationally and assigned
        d1 <= IN_W'(d1_sat);
        d2 <= IN_W'(d2_sat);

        // clock 26-49 is NCDF

        // clock 50
        // x * N(d1)
        // e^-rT * N(d2)
        // e^-rT * N(-d2)
        // x * N(-d1)
        p_xn <= x_d[DELAY_X-1] * nd1;
        p_dn <= dc_d[DELAY_DC-1] * nd2;
        p_pn <= dc_d[DELAY_DC-1] * cd2;
        p_xp <= x_d[DELAY_X-1] * cd1;

        // clock 51
        // call = p_xn - p_dn
        // put = p_pn - p_xp, narrowed
        call_n <= CN_W'(call_s);
        put_n <= PN_W'(put_s);

        // clock 52
        // *K
        p_c <= signed'({1'b0, k_d[DELAY_K-1]}) * call_n;
        p_p <= signed'({1'b0, k_d[DELAY_K-1]}) * put_n;

        // clock 53
        // floor negative numerical noise to 0
        call <= call_raw[PC_W-1] ? '0 : OUT_W'(call_raw);
        put <= put_raw[PP_W-1] ? '0 : OUT_W'(put_raw);

        // delay lines
        s_d[0] <= s;
        s_d[1:DELAY_S-1] <= s_d[0:DELAY_S-2];

        sg_d[0] <= sigma;
        sg_d[1:DELAY_SG-1] <= sg_d[0:DELAY_SG-2];

        tau_d[0] <= tau;
        tau_d[1:DELAY_TAU-1] <= tau_d[0:DELAY_TAU-2];

        r_d[0] <= r;
        for (int i = 1; i < DELAY_R; i++) r_d[i] <= r_d[i-1];

        ht_d[0] <= hv_t;
        ht_d[1:DELAY_HT-1] <= ht_d[0:DELAY_HT-2];

        v2d_d[0] <= v2d;
        v2d_d[1:DELAY_V2D-1] <= v2d_d[0:DELAY_V2D-2];

        x_d[0] <= x;
        x_d[1:DELAY_X-1] <= x_d[0:DELAY_X-2];

        dc_d[0] <= disc;
        dc_d[1:DELAY_DC-1] <= dc_d[0:DELAY_DC-2];

        k_d[0] <= k;
        k_d[1:DELAY_K-1] <= k_d[0:DELAY_K-2];

        if (rst) begin
            vld_u <= '0;
            vld_x <= 1'b0;
            vld_v <= 1'b0;
            vld_num <= 1'b0;
            vld_pd <= 1'b0;
            vld_dw <= 1'b0;
            vld_d <= 1'b0;
            vld_back <= '0;
        end else begin
            vld_u <= {vld_u[0], in_valid};
            vld_x <= ik_valid;
            vld_v <= st_valid;
            vld_num <= lnx_valid;
            vld_pd <= vld_num;
            vld_dw <= vld_pd;
            vld_d <= vld_dw;
            vld_back <= {vld_back[2:0], nd1_valid};
        end
    end

    assign out_valid = vld_back[3];

endmodule
