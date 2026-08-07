// decompose input x = m * 2^e where m in [1,2)
// x == 0 has no decomposition. zero is asserted and m, e are driven to 0.

module lzc_norm #(
    parameter int W = 32,
    parameter int LATENCY = 2,
    localparam int E_W = $clog2(W) + 1
)(
    input logic clk,
    input logic rst,
    input logic in_valid,
    input logic [W-1:0] x,
    output logic out_valid,
    output logic [W-1:0] m,
    output logic signed [E_W-1:0] e,
    output logic zero
);

    localparam int LZ_W = $clog2(W + 1);

    if (LATENCY != 2) begin : g_latency_check
        $error("LATENCY=%0d but the pipeline is 2 deep", LATENCY);
    end

    logic [LZ_W-1:0] lz, lz_q;
    logic [W-1:0] x_q;
    logic [LATENCY-1:0] v;

    // Priority encoder that sets lz = W - 1 - highest index that has 1
    // lz == W means there was none.
    always_comb begin
        lz = LZ_W'(W);
        for (int i = 0; i < W; i++)
            if (x[i]) lz = LZ_W'(W - 1 - i);
    end

    always_ff @(posedge clk) begin
        // clock 0
        lz_q <= lz;
        x_q <= x;

        // clock 1
        m <= x_q << lz_q;
        zero <= (lz_q == LZ_W'(W));
        e <= (lz_q == LZ_W'(W)) ? '0 : E_W'(W - 1 - lz_q);

        if (rst) v <= '0;
        else v <= {v[LATENCY-2:0], in_valid};
    end

    assign out_valid = v[LATENCY-1];

endmodule
