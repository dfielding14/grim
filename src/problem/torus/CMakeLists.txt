# Dump folder
set(OUTPUT_DIR "\"/Users/francoisfoucart/Work/grim/build\"")
set(RESTART "OFF")
set(RESTART_FILE "\"restartfile.h5\"")
set(DUMP_FILE_PREFIX "\"data\"")
set(RESIDUALS_DUMP_FILE_PREFIX "\"residual\"")
set(PROBLEM_DATA "\"${PROBLEM}/problemdata.h\"")

# Time stepping options: EXPLICIT, IMEX or IMPLICIT
set(TIME_STEPPING "IMEX")
set(DT "0.001")
set(DT_DUMP "1.")
set(START_TIME "0.")
set(FINAL_TIME "2000.")
set(START_DUMP_COUNTER "0")
set(COURANT "0.9")
set(MAX_DT_INCREMENT "1.3")

# Domain size. If the problem is 1D then N2 is ignored.
set(COMPUTE_DIM "2")
set(N1 "280")
set(N2 "256")

# The entire global domain is divided into tiles that are small enough to fit
# into the cache of the compute node or an accelerator. This technique
# optimizes cache usage and makes explicit use of the deep memory hierarchies
# which are prevalent in all kinds of processors and accelerators. 
# 
# Caution: 1) The tile sizes need to divide the global domain size in each
#             direction exactly!
#          2) We want a tile size that is as large as the cache. Too large a
#             size and the code can fail silently because local memory in the
#             OpenCL specification is a different amount for different 
#             machines.
#             NEED TO PUT A CHECK FOR THIS.
set(TILE_SIZE_X1 "14")
set(TILE_SIZE_X2 "16")

set(USE_OPENMP   "OFF")

# Physics variables
set(ADIABATIC_INDEX "4./3")
set(FAKE_EMHD "OFF")
set(CONDUCTION "ON")
set(VISCOSITY  "ON")
set(HIGHORDERTERMS_CONDUCTION "ON")
set(HIGHORDERTERMS_VISCOSITY "ON")
set(BACKREACTION "OFF")

set(VISCOSITY_CLOSURE_COEFF "1.")
set(VISCOSITY_ALPHA "1.")

set(CONDUCTION_CLOSURE_COEFF "1.")
set(CONDUCTION_ALPHA "1.")

# Geometry
set(METRIC  "KERRSCHILD")
set(EPS     "1e-5")
set(M       "1.")
set(BH_SPIN "0.9375")
set(H_SLOPE "0.3") # Refinement in the theta direction to better resolve the
                   # disk. H_SLOPE = 1 is no refinement and H_SLOPE -> 0 pushes
                   # all points towards the midplane.

# Domain
set(R_A     "0.98*(1. + sqrt(1. - ${BH_SPIN}*${BH_SPIN}))")
set(R_B     "63.")

set(X1_A  "log(R_A)")
set(X1_B  "log(R_B)")
set(X2_A  "1e-20")
set(X2_B  "1. - 1e-20")

# Boundary conditions
set(PHYSICAL_BOUNDARY_LEFT_EDGE   "OUTFLOW")
set(PHYSICAL_BOUNDARY_RIGHT_EDGE  "OUTFLOW")
set(PHYSICAL_BOUNDARY_TOP_EDGE    "MIRROR")
set(PHYSICAL_BOUNDARY_BOTTOM_EDGE "MIRROR")

# Reconstrution options
# MP5, MONOTONIZED_CENTRAL or MIN_MOD
set(RECONSTRUCTION "MONOTONIZED_CENTRAL")

# Initial condition parameters 
# Everything in units of M
set(R_PRESSURE_MAX          "12.")
set(R_INNER_EDGE            "6.")
set(ADIABAT                 "1e-3")
set(PLASMA_BETA             "1e2")
set(PERTURBATIONS_AMPLITUDE "4e-2")
set(NB_MAGNETIC_LOOPS       "1")

# Floor values
set(RHO_FLOOR         "1e-3")
set(UU_FLOOR          "1e-5")
set(RHO_FLOOR_FALLOFF "-1.5")
set(UU_FLOOR_FALLOFF  "-2.5")
set(RHO_FLOOR_MIN     "1e-15")
set(UU_FLOOR_MIN      "1e-15")
set(GAMMA_MAX         "10.")
set(ADD_WIND_SOURCE   "ON")

message("")
message("#########################")
message("# Configuration options #")
message("#########################")
message("")
message("Problem          : " ${PROBLEM})
message("Output dir       : " ${OUTPUT_DIR})
message("Time stepping    : " ${TIME_STEPPING})
message("Dimensions       : " ${COMPUTE_DIM})
message("Resolution       : " ${N1} " x " ${N2})
message("Tile size        : " ${TILE_SIZE_X1} " x " ${TILE_SIZE_X2})
message("Metric           : " ${METRIC})

message("")
message("#########################")
message("# Geometry options #")
message("#########################")
message("")
message("Metric                   : " ${METRIC})
message("Black hole mass M        : " ${M})
message("Black hole spin BH_SPIN  : " ${BH_SPIN})

message("")
message("##################################")
message("# Domain and boundary conditions #")
message("##################################")
message("")
message("Inner radius                   : " ${R_A})
message("Outer radius                   : " ${R_B})
message("Inner radial boundary          : " ${PHYSICAL_BOUNDARY_LEFT_EDGE})
message("Outer radial boundary          : " ${PHYSICAL_BOUNDARY_RIGHT_EDGE})
message("Bottom theta boundary          : " ${PHYSICAL_BOUNDARY_BOTTOM_EDGE})
message("Top theta boundary             : " ${PHYSICAL_BOUNDARY_TOP_EDGE})
message("Grid refinement factor H_SLOPE : " ${H_SLOPE})

message("")
message("###################")
message("# Physics options #")
message("###################")
message("")
message("Adiabatic index  : " ${ADIABATIC_INDEX})
message("Conduction       : " ${CONDUCTION})
message("Viscosity        : " ${VISCOSITY})
message("Viscosity alpha  : " ${VISCOSITY_ALPHA})
message("Viscosity closure: " ${VISCOSITY_CLOSURE_COEFF})

message("")
message("##################")
message("# Reconstruction #")
message("##################")
message("")
message("Reconstruction  : " ${RECONSTRUCTION})

message("")
message("######################")
message("# Initial conditions #")
message("######################")
message("")
message("R_PRESSURE_MAX   : " ${R_PRESSURE_MAX})
message("R_INNER_EDGE     : " ${R_INNER_EDGE})
message("ADIABAT          : " ${ADIABAT})
message("PLASMA_BETA      : " ${PLASMA_BETA})

message("")
message("##########")
message("# Floors #")
message("##########")
message("")
message("Density floor          : ${RHO_FLOOR} x r^${RHO_FLOOR_FALLOFF}")
message("Internal energy floor  : ${UU_FLOOR} x r^${UU_FLOOR_FALLOFF}")
message("Maximum Lorentz factor : "   ${GAMMA_MAX})
