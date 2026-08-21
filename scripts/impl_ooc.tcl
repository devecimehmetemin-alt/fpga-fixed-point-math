# Out-of-context synth + place + route on the KR260's XCK26, for real post-route timing.
set top   [lindex $argv 0]
set period [lindex $argv 1]
set part  xck26-sfvc784-2LV-c
set root [file normalize [file join [file dirname [info script]] ..]]
set out   [lindex $argv 2]
if {$out eq ""} { set out $root/build }

read_verilog -sv [glob $root/rtl/*.sv]
set toks [lrange $argv 3 end]
set args [list -top $top -part $part -mode out_of_context]
for {set i 0} {$i < [llength $toks]} {incr i} {
    set t [lindex $toks $i]
    if {![string match *=* $t]} { incr i; set t "$t=[lindex $toks $i]" }
    lappend args -generic $t
}
puts "GENERICS: [lrange $args 6 end]"
eval synth_design $args

# the clock port is "clk" on the library modules and "s_axi_aclk" on the wrapper
set clkport [get_ports -quiet clk]
if {$clkport eq ""} { set clkport [get_ports -quiet s_axi_aclk] }
create_clock -name clk -period $period $clkport
set clkname [get_property NAME $clkport]
set io [expr {$period * 0.25}]
set_input_delay -clock clk $io [get_ports -filter "DIRECTION == IN && NAME != $clkname"]
set_output_delay -clock clk $io [get_ports -filter {DIRECTION == OUT}]

opt_design
place_design
phys_opt_design
route_design

file mkdir $out
report_utilization -file $out/${top}_route_util.rpt
report_timing_summary -max_paths 10 -file $out/${top}_route_timing.rpt

set wns [get_property SLACK [get_timing_paths -delay_type max]]
puts "\n======== $top on $part POST-ROUTE, target $period ns ========"
puts "DSP48E2  [llength [get_cells -hier -filter {PRIMITIVE_TYPE =~ *DSP*}]]"
puts "LUT      [llength [get_cells -hier -filter {PRIMITIVE_GROUP == LUT}]]"
puts "FF       [llength [get_cells -hier -filter {PRIMITIVE_GROUP == FLOP_LATCH}]]"
puts "WNS      [format %.3f $wns] ns"
puts "Fmax     [format %.1f [expr {1000.0 / ($period - $wns)}]] MHz"
