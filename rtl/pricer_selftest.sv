module pricer_selftest #(
    parameter int N = 256,
    parameter int IDX_W = 8,
    parameter int SK_W = 32,
    parameter int R_W = 24,
    parameter int SG_W = 24,
    parameter int T_W = 32,
    parameter int OUT_W = 40,
    parameter IN_MEM = "vectors_in.mem",
    parameter EXP_MEM = "vectors_exp.mem"
)(
    input logic clk,
    input logic rst,
    input logic start,
    output logic busy,
    output logic done,
    output logic [IDX_W:0] mismatches,
    output logic [IDX_W-1:0] first_bad,
    output logic [OUT_W-1:0] first_bad_call,
    output logic [OUT_W-1:0] first_bad_put,
    output logic [IDX_W:0] seen
);

    localparam int IN_W = SK_W + SK_W + R_W + SG_W + T_W;
    localparam int EXP_W = 2 * OUT_W;

    logic [IN_W-1:0] rom_in [0:N-1];
    logic [EXP_W-1:0] rom_exp [0:N-1];

    initial begin
        $readmemh(IN_MEM, rom_in);
        $readmemh(EXP_MEM, rom_exp);
    end

    logic [IDX_W-1:0] feed_idx, chk_idx;
    logic feeding, feed_valid;
    logic [IN_W-1:0] in_word;
    logic [EXP_W-1:0] exp_word;

    logic out_valid;
    logic signed [OUT_W-1:0] call, put;

    assign exp_word = rom_exp[chk_idx];
    assign busy = feeding | (seen != (IDX_W+1)'(N));

    bs_pricer_top u_pricer (
        .clk (clk),
        .rst (rst),
        .in_valid (feed_valid),
        .s (in_word[IN_W-1 -: SK_W]),
        .k (in_word[IN_W-SK_W-1 -: SK_W]),
        .r (in_word[T_W+SG_W+R_W-1 -: R_W]),
        .sigma (in_word[T_W+SG_W-1 -: SG_W]),
        .tau (in_word[T_W-1 -: T_W]),
        .out_valid (out_valid),
        .call (call),
        .put (put)
    );

    always_ff @(posedge clk) begin
        feed_valid <= 1'b0;

        if (feeding) begin
            in_word <= rom_in[feed_idx];
            feed_valid <= 1'b1;
            feed_idx <= feed_idx + 1'b1;
            if (feed_idx == IDX_W'(N-1)) feeding <= 1'b0;
        end else if (start & ~done) begin
            feeding <= 1'b1;
            feed_idx <= '0;
        end

        if (out_valid) begin
            seen <= seen + 1'b1;
            chk_idx <= chk_idx + 1'b1;

            if ({$unsigned(call), $unsigned(put)} != exp_word) begin
                if (mismatches == '0) begin
                    first_bad <= chk_idx;
                    first_bad_call <= $unsigned(call);
                    first_bad_put <= $unsigned(put);
                end
                mismatches <= mismatches + 1'b1;
            end

            if (chk_idx == IDX_W'(N-1)) done <= 1'b1;
        end

        if (rst) begin
            feeding <= 1'b0;
            feed_valid <= 1'b0;
            feed_idx <= '0;
            chk_idx <= '0;
            seen <= '0;
            mismatches <= '0;
            first_bad <= '0;
            first_bad_call <= '0;
            first_bad_put <= '0;
            done <= 1'b0;
        end
    end

endmodule
