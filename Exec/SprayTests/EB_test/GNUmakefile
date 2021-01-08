# AMReX
DIM = 3
COMP = llvm
PRECISION = DOUBLE

# Profiling
PROFILE = FALSE
TINY_PROFILE = FALSE
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

#DEFINES += -DCONVERGENCE

# PeleC
USE_REACT = FALSE
USE_EB = TRUE
Eos_dir := Fuego
Reactions_dir := Fuego
Chemistry_Model := heptane_3sp
Transport_dir := Simple

# PeleC-MP
USE_PARTICLES = TRUE
SPRAY_FUEL_NUM = 1

# GNU Make
Bpack := ./Make.package
Blocs := .
include $(PELEC_HOME)/ExecCpp/Make.PeleC
