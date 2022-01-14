#
# Build configuration
#

DIM             = 2
COMP            = gnu
DEBUG           = FALSE
USE_MPI         = FALSE
USE_OMP         = FALSE
PRECISION       = DOUBLE

BL_NO_FORT = TRUE
#WARN_ALL = TRUE

USE_CUDA        = FALSE

TINY_PROFILE    = TRUE

USE_EB = FALSE

VERBOSE = FALSE

USE_SUNDIALS = TRUE
Chemistry_Model := decane_3sp

# This sets the EOS directory in $(PELE_PHYSICS_HOME)/Eos
Eos_Model     := Fuego

# This sets the transport directory in $(PELE_PHYSICS_HOME)/Transport
Transport_Model := Simple

USE_PARTICLES = TRUE
SPRAY_FUEL_NUM = 1
CEXE_sources += SprayParticlesInitInsert.cpp

Blocs   := .

include $(PELELM_HOME)/Tools/Make/Make.PeleLM
