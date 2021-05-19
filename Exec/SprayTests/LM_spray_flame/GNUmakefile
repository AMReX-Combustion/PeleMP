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
DEBUG           = TRUE
PRECISION       = DOUBLE
VERBOSE         = FALSE
TINY_PROFILE    = FALSE

# CVODE
USE_SUNDIALS_PP = TRUE
USE_KLU_PP      = FALSE

# PeleLM options
DO_2S_CONVERGENCE=FALSE

Chemistry_Model := decane_3sp

# This sets the EOS directory in $(PELE_PHYSICS_HOME)/Eos
Eos_dir     := Fuego

# This sets the network directory in $(PELE_PHYSICS_HOME)/Reactions
Reactions_dir := Fuego

# This sets the transport directory in $(PELE_PHYSICS_HOME)/Transport
Transport_dir := Simple

USE_PARTICLES = TRUE
SPRAY_FUEL_NUM = 1
CEXE_sources += SprayParticlesInitInsert.cpp

Blocs   := .

CEXE_sources +=
F90EXE_sources +=
CEXE_headers += 
FEXE_headers += 

include $(PELELM_HOME)/Tools/Make/Make.PeleLM

ifeq ($(USE_CUDA),TRUE)
CXXFLAGS += -Xptxas --disable-optimizer-constants
endif

