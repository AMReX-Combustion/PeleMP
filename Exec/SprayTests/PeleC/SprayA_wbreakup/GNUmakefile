# AMReX
DIM = 3
COMP = gnu
PRECISION = DOUBLE

# Profiling
PROFILE = FALSE
TINY_PROFILE = TRUE
COMM_PROFILE = FALSE
TRACE_PROFILE = FALSE
MEM_PROFILE = FALSE
USE_GPROF = FALSE

# Performance
USE_MPI = TRUE
USE_OMP = FALSE
USE_CUDA = FALSE
USE_HIP = FALSE
PELE_USE_MAGMA = TRUE

# Debugging
DEBUG = FALSE
FSANITIZER = FALSE
THREAD_SANITIZER = FALSE

# PeleC
Eos_Model := Fuego
Chemistry_Model := dodecane_lu
Transport_Model := Simple

# PeleMP
USE_PARTICLES = TRUE
SPRAY_FUEL_NUM = 1

# GNU Make
Bpack := ./Make.package
Blocs := .
include $(PELEC_HOME)/Exec/Make.PeleC
