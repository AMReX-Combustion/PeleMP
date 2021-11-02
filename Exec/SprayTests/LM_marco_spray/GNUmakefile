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
USE_CUDA        = TRUE
USE_HIP         = FALSE
USE_FORTRAN_INTERFACE = FALSE

# MISC options
DEBUG           = FALSE
PRECISION       = DOUBLE
VERBOSE         = FALSE
TINY_PROFILE    = TRUE

# Ascent
USE_ASCENT      = FALSE
USE_CONDUIT     = FALSE

# CVODE
#USE_SUNDIALS_PP = TRUE
USE_KLU_PP      = FALSE

# PeleLM options
DO_2S_CONVERGENCE=FALSE

Chemistry_Model := dodecane_lu
###Chemistry_Model := decane_3sp

# INCLUDE_LOCATIONS +=${PELE_PHYSICS_HOME}/ThirdParty/INSTALL/gcc.CUDA/include/
# LIBRARY_LOCATIONS +=${PELE_PHYSICS_HOME}/ThirdParty/INSTALL/gcc.CUDA/lib/
#
# This sets the EOS directory in $(PELE_PHYSICS_HOME)/Eos
Eos_Model   := Fuego

# This sets the network directory in $(PELE_PHYSICS_HOME)/Reactions

# This sets the transport directory in $(PELE_PHYSICS_HOME)/Transport
Transport_Model := Simple

USE_PARTICLES = TRUE
SPRAY_FUEL_NUM = 1
CEXE_sources += SprayParticlesInitInsert.cpp

Blocs   := .

# ifeq ($(USE_ASCENT), TRUE)

# DEFINES += -DASCENT
# DEFINES += -DUSE_ASCENT
# #ASCENT_INSTALL_DIR=/gpfs/alpine/world-shared/csc340/cokurt/ascent_gpu_9_21/ascent/install
# #/gpfs/alpine/world-shared/csc340/cokurt/ascent_cokurt_cuda/install/
# ASCENT_INSTALL_DIR=/gpfs/alpine/world-shared/csc340/cokurt/ascent_genten/ascent/install/

# include $(ASCENT_INSTALL_DIR)/share/ascent/ascent_config.mk

# include $(AMREX_HOME)/Src/Extern/Conduit/Make.package
# ASCENT_LINK_FLAGS = $(subst -pthread,, $(ASCENT_MPI_LIB_FLAGS))
# INCLUDE_LOCATIONS += $(ASCENT_INCLUDE_FLAGS)
# INCLUDE_LOCATIONS += $(AMREX_HOME)/Src/Extern/Conduit
# INCLUDE_LOCATIONS += ${ASCENT_CONDUIT_DIR}/include/conduit/
# LIBRARIES += $(ASCENT_LINK_FLAGS)

# endif

# CEXE_sources +=
# F90EXE_sources +=
# CEXE_headers += 
# FEXE_headers += 

include $(PELELM_HOME)/Tools/Make/Make.PeleLM

ifeq ($(USE_CUDA),TRUE)
CXXFLAGS += -Xptxas --disable-optimizer-constants
endif

