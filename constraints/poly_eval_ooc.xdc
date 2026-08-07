# Timing constraints for out-of-context characterisation of poly_eval.
#
# This is not a board constraint file. There are no pin assignments, because
# poly_eval is never a top-level design -- it exists to be instantiated. The
# only job here is to give report_timing_summary a clock, so that Fmax is a
# measurement rather than "inf".
#
#   read_verilog -sv rtl/poly_eval.sv
#   synth_design -top poly_eval -part xc7a200tsbg484-1 -mode out_of_context
#   read_xdc constraints/poly_eval_ooc.xdc
#   opt_design
#   report_timing_summary -max_paths 5

# Sweep this to find Fmax: tighten until WNS goes negative, then back off to
# the smallest period that still closes. 3.000 ns is 333 MHz and is meant to
# fail on the first run -- a constraint that passes immediately has told you
# nothing about where the critical path is.
set PERIOD 3.000

create_clock -name clk -period $PERIOD [get_ports clk]

# Budget a quarter period either side for the logic that will eventually
# surround this module. Without these, the I/O paths are skipped entirely and
# the summary reports zero failing endpoints while having checked almost
# nothing -- the same silent pass that gave WNS = inf.
set IO_BUDGET [expr {$PERIOD * 0.25}]

set_input_delay  -clock clk $IO_BUDGET \
    [remove_from_collection [all_inputs] [get_ports clk]]
set_output_delay -clock clk $IO_BUDGET [all_outputs]

# Keep the boundary registers out of the IOBs. In out-of-context mode there
# are no IO buffers and these lines do nothing, but on a top-level run Vivado
# will happily pack x_r and p into IOB flip-flops -- which blocks synthesis
# from absorbing them into the DSP48 A/B and P registers, inflating the fabric
# FF count and moving the critical path somewhere it will not be in the real
# design.
set_property IOB FALSE [all_inputs]
set_property IOB FALSE [all_outputs]
