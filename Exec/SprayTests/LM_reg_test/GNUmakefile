#
# Build configuration
#

DIM             = 2
COMP            = gnu
DEBUG           = FALSE
USE_MPI         = FALSE
USE_OMP         = FALSE
PRECISION       = DOUBLE

#BL_NO_FORT = TRUE
#WARN_ALL = TRUE

USE_CUDA        = FALSE

TINY_PROFILE    = TRUE

USE_EB = FALSE

VERBOSE = FALSE
USE_SUNDIALS_PP = TRUE
USE_KLU_PP = FALSE

# Enable ht.velocity_plotfile option
#USE_VELOCITY    = FALSE

Chemistry_Model := decane_3sp

# This sets the EOS directory in $(PELE_PHYSICS_HOME)/Eos
Eos_dir     := Fuego

# This sets the network directory in $(PELE_PHYSICS_HOME)/Reactions
Reactor_dir := cvode

# This sets the transport directory in $(PELE_PHYSICS_HOME)/Transport
Transport_dir := Simple

USE_PARTICLES = TRUE
SPRAY_FUEL_NUM = 1
CEXE_sources += SprayParticlesInitInsert.cpp

Blocs   := .

include $(PELELM_HOME)/Tools/Make/Make.PeleLM
