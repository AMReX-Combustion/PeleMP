# AMReX
DIM = 2
COMP = gnu
PRECISION = DOUBLE

BL_NO_FORT = TRUE

# Profiling
PROFILE = FALSE
TINY_PROFILE = TRUE
COMM_PROFILE = FALSE
TRACE_PROFILE = FALSE
MEM_PROFILE = FALSE
USE_GPROF = FALSE

# Performance
USE_MPI = FALSE
USE_OMP = FALSE
USE_CUDA = FALSE
USE_HIP = FALSE

# Debugging
DEBUG = FALSE
FSANITIZER = FALSE
THREAD_SANITIZER = FALSE

# PeleC
USE_REACT = FALSE
USE_EB = FALSE
Eos_dir := GammaLaw
Reactions_dir := Null
Transport_dir := Constant

# PeleC-MP
USE_PARTICLES = TRUE
SPRAY_FUEL_NUM = 1
# Set to true if we want to use StructOfArray particles data
#USE_SPRAY_SOA = TRUE

# GNU Make
Bpack := ./Make.package
Blocs := .
include $(PELEC_HOME)/ExecCpp/Make.PeleC
