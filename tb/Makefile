SIM ?= verilator
TOPLEVEL_LANG ?= verilog

VERILOG_SOURCES += $(PWD)/../rtl/poly_eval.sv
VERILOG_SOURCES += $(PWD)/../rtl/exp_unit.sv
TOPLEVEL = exp_unit
COCOTB_TEST_MODULES = tb_exp_unit
DUT ?= exp_unit

# Keep the C++ compile off /mnt/c -- it is roughly an order of magnitude slower.
SIM_BUILD = $(HOME)/mc-build/$(DUT)

EXTRA_ARGS += --timing -Wall

# make WAVES=1 -- opt-in, since tracing slows both the build and the run.
# --trace-structs keeps the packed coefficient array readable as a[0]..a[5].
ifeq ($(WAVES),1)
EXTRA_ARGS += --trace --trace-structs
endif

# Override on the command line to sweep configurations, e.g.
#   make X_W=18 X_F=17 C_W=25 C_F=23 W_W=24 W_F=22
X_W ?= 18
X_F ?= 17
L_W ?= 18
L_F ?= 16
T_W ?= 20
T_F ?= 18
C_W ?= 18
C_F ?= 16
W_W ?= 18
W_F ?= 16
Y_W ?= 32
Y_F ?= 31
MAX_SHIFT ?= 1
LATENCY ?= 12

COMPILE_ARGS += -GX_W=$(X_W) -GX_F=$(X_F) -GC_W=$(C_W) -GC_F=$(C_F) -GW_W=$(W_W) -GW_F=$(W_F)
COMPILE_ARGS += -GL_W=$(L_W) -GL_F=$(L_F) -GT_W=$(T_W) -GT_F=$(T_F) -GY_W=$(Y_W) -GY_F=$(Y_F)
COMPILE_ARGS += -GMAX_SHIFT=$(MAX_SHIFT) -GLATENCY=$(LATENCY)

# Same numbers to the testbench, so the formats are defined in one place only.
export X_W X_F L_W L_F T_W T_F C_W C_F W_W W_F Y_W Y_F MAX_SHIFT LATENCY

export PYTHONPATH := $(PWD)/../model:$(PYTHONPATH)

include $(shell cocotb-config --makefiles)/Makefile.sim
