# Out-of-context synthesis of one module: DSP/LUT/FF counts and Fmax.
#
# No board, no pins, no top level -- the point is to characterise a block in
# isolation. Tighten the period until WNS goes negative, then back off: a
# constraint that passes first try tells you nothing about the critical path.
#
#   vivado -mode batch -source scripts/synth_ooc.tcl -tclargs poly_eval
#   vivado -mode batch -source scripts/synth_ooc.tcl -tclargs exp_unit 3.0
#   vivado -mode batch -source scripts/synth_ooc.tcl -tclargs exp_unit 3.0 X_W=24 X_F=19

set top [lindex $argv 0]
if {$top eq ""} { error "usage: -tclargs <top> \[period_ns\] \[NAME=VALUE ...\]" }

set period 3.0
if {[llength $argv] > 1} { set period [lindex $argv 1] }

set part xc7a200tsbg484-1
set root [file normalize [file join [file dirname [info script]] ..]]

read_verilog -sv [glob $root/rtl/*.sv]

# -tclargs splits NAME=VALUE into two tokens, which silently binds every
# parameter to 0 and collapses the design. Rejoin them.
set toks [lrange $argv 2 end]
set args [list -top $top -part $part -mode out_of_context]
for {set i 0} {$i < [llength $toks]} {incr i} {
    set t [lindex $toks $i]
    if {![string match *=* $t]} {
        incr i
        set t "$t=[lindex $toks $i]"
    }
    lappend args -generic $t
}
eval synth_design $args

# Constrain after synthesis: WNS = inf means no clock was defined, not that
# timing is met. A quarter period each side keeps the I/O paths from being
# skipped entirely.
create_clock -name clk -period $period [get_ports clk]
set io [expr {$period * 0.25}]
set_input_delay  -clock clk $io [get_ports -filter {DIRECTION == IN && NAME != clk}]
set_output_delay -clock clk $io [get_ports -filter {DIRECTION == OUT}]

opt_design

file mkdir $root/build
report_utilization -hierarchical -file $root/build/${top}_util.rpt
report_timing_summary -max_paths 5 -file $root/build/${top}_timing.rpt

set wns [get_property SLACK [get_timing_paths -delay_type max]]
puts "\n======== $top on $part, target [format %.3f $period] ns ========"
puts "DSP48E1  [llength [get_cells -hier -filter {PRIMITIVE_TYPE =~ *DSP*}]]"
puts "LUT      [llength [get_cells -hier -filter {PRIMITIVE_GROUP == LUT}]]"
puts "FF       [llength [get_cells -hier -filter {PRIMITIVE_GROUP == FLOP_LATCH}]]"
puts "SRL      [llength [get_cells -hier -filter {PRIMITIVE_SUBGROUP == SRL}]]"
puts "WNS      [format %.3f $wns] ns"
puts "Fmax     [format %.1f [expr {1000.0 / ($period - $wns)}]] MHz"
puts "reports  build/${top}_util.rpt  build/${top}_timing.rpt"
