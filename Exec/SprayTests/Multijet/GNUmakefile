# AMReX
DIM = 3
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
USE_REACT = TRUE
USE_EB = FALSE
Eos_Model := Fuego
Reactor_dir := cvode
Chemistry_Model := dodecane_lu
Transport_Model := Simple

# PeleC-MP
USE_PARTICLES = TRUE
SPRAY_FUEL_NUM = 1
#USE_SPRAY_SOA = TRUE

# GNU Make
Bpack := ./Make.package
Blocs := .
include $(PELEC_HOME)/Exec/Make.PeleC
