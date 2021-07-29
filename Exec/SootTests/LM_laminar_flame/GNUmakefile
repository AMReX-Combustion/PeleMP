# AMReX
DIM = 2
COMP = llvm
PRECISION = DOUBLE

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
USE_SUNDIALS_PP = TRUE
USE_EB = FALSE
Eos_dir := Fuego
Reactor_dir := cvode
Chemistry_Model := SootReaction
Transport_dir := Simple

# PeleC-MP
USE_SOOT = TRUE
# If this is changed, must run a make clean and rerun make
NUM_SOOT_MOMENTS = 3
#NUM_SOOT_MOMENTS = 6

# GNU Make
Bpack := ./Make.package
Blocs := .
include $(PELELM_HOME)/Tools/Make/Make.PeleLM
