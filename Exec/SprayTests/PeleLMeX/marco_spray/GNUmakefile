

# AMReX
DIM             = 3
DEBUG           = FALSE
PRECISION       = DOUBLE
VERBOSE         = FALSE
TINY_PROFILE    = FALSE

# Compilation
COMP            = gnu
USE_MPI         = TRUE
USE_OMP         = FALSE
USE_CUDA        = TRUE
USE_HIP         = FALSE

# PeleLMeX
USE_EFIELD      = FALSE

# Uncomment for Blaze
#DEFINES+=-DAMREX_GPUS_PER_NODE=2

# PelePhysics
Chemistry_Model = dodecane_lu
Eos_Model	= Fuego
Transport_Model = Simple

USE_PARTICLES = TRUE
SPRAY_FUEL_NUM = 1

include $(PELELMEX_HOME)/Utils/Make.PeleLMeX
