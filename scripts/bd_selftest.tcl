set root [file normalize [file join [file dirname [info script]] ..]]
set proj $root/build/kr260_selftest
set part xck26-sfvc784-2LV-c

create_project kr260_selftest $proj -part $part -force
set_property board_part [lindex [lsort [get_board_parts -quiet *kr260_som*]] end] [current_project]

add_files -norecurse [glob $root/rtl/*.sv]
add_files -norecurse [list $root/tb/vectors_in.mem $root/tb/vectors_exp.mem]
set_property file_type SystemVerilog [get_files *.sv]
update_compile_order -fileset sources_1

create_bd_design "top"

set ps [create_bd_cell -type ip -vlnv xilinx.com:ip:zynq_ultra_ps_e ps]
apply_bd_automation -rule xilinx.com:bd_rule:zynq_ultra_ps_e -config {apply_board_preset "1"} $ps
set_property -dict [list \
    CONFIG.PSU__USE__M_AXI_GP0 {0} \
    CONFIG.PSU__USE__M_AXI_GP1 {0} \
    CONFIG.PSU__USE__M_AXI_GP2 {0} \
    CONFIG.PSU__FPGA_PL0_ENABLE {1} \
    CONFIG.PSU__CRL_APB__PL0_REF_CTRL__FREQMHZ {100} \
] $ps

set dut [create_bd_cell -type module -reference pricer_selftest st]

set vio [create_bd_cell -type ip -vlnv xilinx.com:ip:vio vio]
set_property -dict [list \
    CONFIG.C_NUM_PROBE_IN {6} \
    CONFIG.C_NUM_PROBE_OUT {2} \
    CONFIG.C_PROBE_IN0_WIDTH {1} \
    CONFIG.C_PROBE_IN1_WIDTH {9} \
    CONFIG.C_PROBE_IN2_WIDTH {8} \
    CONFIG.C_PROBE_IN3_WIDTH {40} \
    CONFIG.C_PROBE_IN4_WIDTH {40} \
    CONFIG.C_PROBE_IN5_WIDTH {9} \
    CONFIG.C_PROBE_OUT0_WIDTH {1} \
    CONFIG.C_PROBE_OUT1_WIDTH {1} \
    CONFIG.C_PROBE_OUT0_INIT_VAL {0x0} \
    CONFIG.C_PROBE_OUT1_INIT_VAL {0x1} \
] $vio

connect_bd_net [get_bd_pins ps/pl_clk0] [get_bd_pins st/clk] [get_bd_pins vio/clk]
connect_bd_net [get_bd_pins vio/probe_out0] [get_bd_pins st/start]
connect_bd_net [get_bd_pins vio/probe_out1] [get_bd_pins st/rst]
connect_bd_net [get_bd_pins st/done] [get_bd_pins vio/probe_in0]
connect_bd_net [get_bd_pins st/mismatches] [get_bd_pins vio/probe_in1]
connect_bd_net [get_bd_pins st/first_bad] [get_bd_pins vio/probe_in2]
connect_bd_net [get_bd_pins st/first_bad_call] [get_bd_pins vio/probe_in3]
connect_bd_net [get_bd_pins st/first_bad_put] [get_bd_pins vio/probe_in4]
connect_bd_net [get_bd_pins st/seen] [get_bd_pins vio/probe_in5]

validate_bd_design
save_bd_design
make_wrapper -files [get_files top.bd] -top
add_files -norecurse [glob $proj/*.gen/sources_1/bd/top/hdl/top_wrapper.v]
set_property top top_wrapper [current_fileset]

launch_runs impl_1 -to_step write_bitstream -jobs 8
wait_on_run impl_1
puts "BITSTREAM $proj/kr260_selftest.runs/impl_1/top_wrapper.bit"
puts "PROBES    $proj/kr260_selftest.runs/impl_1/top_wrapper.ltx"
