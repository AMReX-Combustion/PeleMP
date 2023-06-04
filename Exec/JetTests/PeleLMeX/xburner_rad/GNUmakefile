

# AMReX
DIM             = 3
DEBUG           = FALSE
PRECISION       = DOUBLE
VERBOSE         = TRUE
TINY_PROFILE    = TRUE

# Compilation
COMP            = clang
USE_MPI         = TRUE
USE_OMP         = FALSE
USE_CUDA        = FALSE
USE_HIP         = TRUE

# PeleLMeX
USE_EFIELD      = FALSE

# PelePhysics
Chemistry_Model = SootReaction
Eos_Model	= Fuego
Transport_Model = Simple

USE_SOOT = TRUE
NUM_SOOT_MOMENTS = 3

# PeleRad
USE_PELERAD = TRUE

include $(PELELMEX_HOME)/Utils/Make.PeleLMeX
