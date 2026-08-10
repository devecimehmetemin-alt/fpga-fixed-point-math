SIM ?= verilator
TOPLEVEL_LANG ?= verilog

# make DUT=lzc_norm  --  selects the module, its testbench and its parameters.
DUT ?= exp_unit
TOPLEVEL = $(DUT)
COCOTB_TEST_MODULES = tb_$(DUT)

VERILOG_SOURCES += $(PWD)/../rtl/$(DUT).sv
ifeq ($(DUT),exp_unit)
VERILOG_SOURCES += $(PWD)/../rtl/poly_eval.sv
endif
ifeq ($(DUT),ln_unit)
VERILOG_SOURCES += $(PWD)/../rtl/lzc_norm.sv $(PWD)/../rtl/poly_eval.sv
endif
ifeq ($(DUT),recip_unit)
VERILOG_SOURCES += $(PWD)/../rtl/lzc_norm.sv
endif
ifeq ($(DUT),sqrt_unit)
VERILOG_SOURCES += $(PWD)/../rtl/lzc_norm.sv
endif

EXTRA_ARGS += --timing -Wall

# make WAVES=1 -- opt-in, since tracing slows both the build and the run.
# --trace-structs keeps the packed coefficient array readable as a[0]..a[5].
ifeq ($(WAVES),1)
EXTRA_ARGS += --trace --trace-structs
endif

# Override any of these on the command line to sweep configurations. The same
# numbers go to Verilator as -G and to the testbench as environment variables,
# so the formats are defined in one place only.

ifeq ($(DUT),poly_eval)
X_W ?= 18
X_F ?= 16
C_W ?= 18
C_F ?= 16
W_W ?= 18
W_F ?= 15
COMPILE_ARGS += -GX_W=$(X_W) -GX_F=$(X_F) -GC_W=$(C_W) -GC_F=$(C_F) -GW_W=$(W_W) -GW_F=$(W_F)
export X_W X_F C_W C_F W_W W_F
endif

ifeq ($(DUT),exp_unit)
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
export X_W X_F L_W L_F T_W T_F C_W C_F W_W W_F Y_W Y_F MAX_SHIFT LATENCY
endif

ifeq ($(DUT),lzc_norm)
W ?= 12
LATENCY ?= 2
COMPILE_ARGS += -GW=$(W) -GLATENCY=$(LATENCY)
export W LATENCY
endif

ifeq ($(DUT),ln_unit)
W ?= 32
F ?= 28
X_W ?= 18
X_F ?= 17
C_W ?= 18
C_F ?= 17
W_W ?= 18
W_F ?= 16
L_W ?= 18
Y_W ?= 24
LATENCY ?= 11
COMPILE_ARGS += -GW=$(W) -GF=$(F) -GX_W=$(X_W) -GX_F=$(X_F) -GC_W=$(C_W) -GC_F=$(C_F)
COMPILE_ARGS += -GW_W=$(W_W) -GW_F=$(W_F) -GL_W=$(L_W) -GY_W=$(Y_W) -GLATENCY=$(LATENCY)
export W F X_W X_F C_W C_F W_W W_F L_W Y_W LATENCY
# m's low W-1-X_F bits are dead by construction: z takes m[W-2 -: X_F] and the
# leading one is implicit. Put a lint_off pragma on the declaration in the RTL
# to drop this.
EXTRA_ARGS += -Wno-UNUSEDSIGNAL
endif

ifeq ($(DUT),recip_unit)
A_W ?= 32
A_F ?= 28
Y_W ?= 24
Y_F ?= 22
IDX_W ?= 9
S_F ?= 17
M_W ?= 25
R_W ?= 18
E_MIN ?= 28
E_MAX ?= 29
LATENCY ?= 10
COMPILE_ARGS += -GA_W=$(A_W) -GA_F=$(A_F) -GY_W=$(Y_W) -GY_F=$(Y_F) -GIDX_W=$(IDX_W)
COMPILE_ARGS += -GS_F=$(S_F) -GM_W=$(M_W) -GR_W=$(R_W) -GE_MIN=$(E_MIN) -GE_MAX=$(E_MAX)
COMPILE_ARGS += -GLATENCY=$(LATENCY)
export A_W A_F Y_W Y_F IDX_W S_F M_W R_W E_MIN E_MAX LATENCY
# m's bits below m_q's slice are dead by construction. Put a lint_off pragma on
# the declaration in the RTL to drop this.
EXTRA_ARGS += -Wno-UNUSEDSIGNAL
endif

ifeq ($(DUT),sqrt_unit)
A_W ?= 32
A_F ?= 28
Y_W ?= 24
Y_F ?= 21
IDX_W ?= 10
S_F ?= 17
M_W ?= 25
R_W ?= 18
E_MIN ?= 28
E_MAX ?= 29
LATENCY ?= 12
COMPILE_ARGS += -GA_W=$(A_W) -GA_F=$(A_F) -GY_W=$(Y_W) -GY_F=$(Y_F) -GIDX_W=$(IDX_W)
COMPILE_ARGS += -GS_F=$(S_F) -GM_W=$(M_W) -GR_W=$(R_W) -GE_MIN=$(E_MIN) -GE_MAX=$(E_MAX)
COMPILE_ARGS += -GLATENCY=$(LATENCY)
export A_W A_F Y_W Y_F IDX_W S_F M_W R_W E_MIN E_MAX LATENCY
# m's bits below m_q's slice, and lzc_norm's zero flag, are dead by
# construction. sqrt(0) falls out of the underflow branch on its own.
EXTRA_ARGS += -Wno-UNUSEDSIGNAL
endif

# Keep the C++ compile off /mnt/c -- it is roughly an order of magnitude slower.
# The directory has to carry the parameters. cocotb rebuilds when a source file
# changes, not when a -G value does, so a shared build directory hands a sweep
# the previous run's netlist while the model uses the new numbers -- the RTL
# looks wrong and is not. Must come after the blocks above, since := reads
# COMPILE_ARGS now.
BUILD_TAG := $(shell echo "$(COMPILE_ARGS) $(EXTRA_ARGS)" | md5sum | cut -c1-8)
SIM_BUILD = $(HOME)/mc-build/$(DUT)-$(BUILD_TAG)

export PYTHONPATH := $(PWD)/../model:$(PYTHONPATH)

include $(shell cocotb-config --makefiles)/Makefile.sim
