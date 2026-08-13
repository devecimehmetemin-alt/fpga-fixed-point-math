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
ifeq ($(DUT),norm_cdf)
VERILOG_SOURCES += $(PWD)/../rtl/lzc_norm.sv $(PWD)/../rtl/poly_eval.sv
VERILOG_SOURCES += $(PWD)/../rtl/recip_unit.sv $(PWD)/../rtl/exp_unit.sv
endif
ifeq ($(DUT),bs_pricer_top)
VERILOG_SOURCES += $(PWD)/../rtl/lzc_norm.sv $(PWD)/../rtl/poly_eval.sv
VERILOG_SOURCES += $(PWD)/../rtl/recip_unit.sv $(PWD)/../rtl/sqrt_unit.sv
VERILOG_SOURCES += $(PWD)/../rtl/ln_unit.sv $(PWD)/../rtl/exp_unit.sv
VERILOG_SOURCES += $(PWD)/../rtl/norm_cdf.sv
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

ifeq ($(DUT),norm_cdf)
IN_W ?= 26
IN_F ?= 22
CLAMP ?= 6.0
P_W ?= 18
P_F ?= 17
T_W ?= 25
T_F ?= 23
C_W ?= 18
C_F ?= 17
W_W ?= 25
W_F ?= 24
U_W ?= 28
U_F ?= 22
PHI_W ?= 32
PHI_F ?= 31
N_W ?= 26
N_F ?= 25
IDX_W ?= 10
EC_W ?= 23
EC_F ?= 21
EW_W ?= 24
EW_F ?= 22
LATENCY ?= 23
COMPILE_ARGS += -GIN_W=$(IN_W) -GIN_F=$(IN_F) -GCLAMP=$(CLAMP) -GP_W=$(P_W) -GP_F=$(P_F)
COMPILE_ARGS += -GT_W=$(T_W) -GT_F=$(T_F) -GC_W=$(C_W) -GC_F=$(C_F) -GW_W=$(W_W) -GW_F=$(W_F)
COMPILE_ARGS += -GU_W=$(U_W) -GU_F=$(U_F) -GPHI_W=$(PHI_W) -GPHI_F=$(PHI_F)
COMPILE_ARGS += -GN_W=$(N_W) -GN_F=$(N_F) -GIDX_W=$(IDX_W)
COMPILE_ARGS += -GEC_W=$(EC_W) -GEC_F=$(EC_F) -GEW_W=$(EW_W) -GEW_F=$(EW_F)
COMPILE_ARGS += -GLATENCY=$(LATENCY)
export IN_W IN_F CLAMP P_W P_F T_W T_F C_W C_F W_W W_F U_W U_F
export PHI_W PHI_F N_W N_F IDX_W EC_W EC_F EW_W EW_F LATENCY
# recip_unit's dead m bits again, plus two outputs that are dead by
# construction here: arg_ovf (arg >= 1, so no zero and no saturation) and
# exp_unit's out_valid (it has to coincide with poly_eval's). Both are
# assertions, not logic.
EXTRA_ARGS += -Wno-UNUSEDSIGNAL
endif

ifeq ($(DUT),bs_pricer_top)
SK_W ?= 32
SK_F ?= 22
T_W ?= 32
T_F ?= 28
SG_W ?= 24
SG_F ?= 22
R_W ?= 24
R_F ?= 22
IK_W ?= 32
IK_F ?= 32
X_W ?= 36
X_F ?= 30
LX_W ?= 28
LX_F ?= 24
ST_W ?= 28
ST_F ?= 25
V_W ?= 28
V_F ?= 25
IV_W ?= 32
IV_F ?= 19
HV_W ?= 26
HV_F ?= 24
HT_W ?= 28
HT_F ?= 24
NM_W ?= 29
NM_F ?= 24
D_W ?= 40
IN_W ?= 26
IN_F ?= 22
CLAMP ?= 6.0
N_W ?= 26
N_F ?= 25
U_W ?= 28
U_F ?= 22
DC_W ?= 32
DC_F ?= 31
CN_W ?= 32
CN_F ?= 26
PN_W ?= 28
PN_F ?= 26
OUT_W ?= 40
OUT_F ?= 22
K_IDX ?= 12
V_IDX ?= 10
EC_W ?= 23
EC_F ?= 21
EW_W ?= 24
EW_F ?= 22
LATENCY ?= 53
COMPILE_ARGS += -GSK_W=$(SK_W) -GSK_F=$(SK_F) -GT_W=$(T_W) -GT_F=$(T_F)
COMPILE_ARGS += -GSG_W=$(SG_W) -GSG_F=$(SG_F) -GR_W=$(R_W) -GR_F=$(R_F)
COMPILE_ARGS += -GIK_W=$(IK_W) -GIK_F=$(IK_F) -GX_W=$(X_W) -GX_F=$(X_F)
COMPILE_ARGS += -GLX_W=$(LX_W) -GLX_F=$(LX_F) -GST_W=$(ST_W) -GST_F=$(ST_F)
COMPILE_ARGS += -GV_W=$(V_W) -GV_F=$(V_F) -GIV_W=$(IV_W) -GIV_F=$(IV_F)
COMPILE_ARGS += -GHV_W=$(HV_W) -GHV_F=$(HV_F) -GHT_W=$(HT_W) -GHT_F=$(HT_F)
COMPILE_ARGS += -GNM_W=$(NM_W) -GNM_F=$(NM_F) -GD_W=$(D_W)
COMPILE_ARGS += -GIN_W=$(IN_W) -GIN_F=$(IN_F) -GCLAMP=$(CLAMP)
COMPILE_ARGS += -GN_W=$(N_W) -GN_F=$(N_F) -GU_W=$(U_W) -GU_F=$(U_F)
COMPILE_ARGS += -GDC_W=$(DC_W) -GDC_F=$(DC_F) -GCN_W=$(CN_W) -GCN_F=$(CN_F)
COMPILE_ARGS += -GPN_W=$(PN_W) -GPN_F=$(PN_F) -GOUT_W=$(OUT_W) -GOUT_F=$(OUT_F)
COMPILE_ARGS += -GK_IDX=$(K_IDX) -GV_IDX=$(V_IDX)
COMPILE_ARGS += -GEC_W=$(EC_W) -GEC_F=$(EC_F) -GEW_W=$(EW_W) -GEW_F=$(EW_F)
COMPILE_ARGS += -GLATENCY=$(LATENCY)
export SK_W SK_F T_W T_F SG_W SG_F R_W R_F IK_W IK_F X_W X_F LX_W LX_F
export ST_W ST_F V_W V_F IV_W IV_F HV_W HV_F HT_W HT_F NM_W NM_F D_W
export IN_W IN_F CLAMP N_W N_F U_W U_F DC_W DC_F CN_W CN_F PN_W PN_F
export OUT_W OUT_F K_IDX V_IDX EC_W EC_F EW_W EW_F LATENCY
# The dead m bits in both recip instances and in sqrt, plus four sub-module
# outputs that are dead by construction here: recip_k's ovf (K >= 10), ln's zero
# (x >= 1/30), and norm_cdf's two clamp-related internals. All assertions
# rather than logic.
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
