#PELELM_HOME = ${PELELM_HOME} #../../..
#SUBMODS = ${PELELM_HOME}/Submodules
#AMREX_HOME         ?= ${SUBMODS}/amrex
#IAMR_HOME          ?= ${SUBMODS}/IAMR
#PELE_PHYSICS_HOME  ?= ${SUBMODS}/PelePhysics

#
# Build configuration
#

# AMREX options
DIM             = 3

# Compiler / parrallel paradigms
COMP            = gnu
USE_MPI         = TRUE
USE_OMP         = FALSE
USE_CUDA        = FALSE
USE_HIP         = FALSE

# MISC options
DEBUG           = FALSE
PRECISION       = DOUBLE
VERBOSE         = FALSE
TINY_PROFILE    = TRUE

Chemistry_Model := decane_3sp

# This sets the EOS directory in $(PELE_PHYSICS_HOME)/Eos
Eos_Model     := Fuego

# This sets the transport directory in $(PELE_PHYSICS_HOME)/Transport
Transport_Model := Simple

USE_PARTICLES = TRUE
SPRAY_FUEL_NUM = 1

include $(PELELM_HOME)/Tools/Make/Make.PeleLM

