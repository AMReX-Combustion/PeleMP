PRECISION  = DOUBLE
PROFILE    = FALSE
DEBUG      = FALSE

DIM        = 3

COMP	   = gcc

USE_MPI    = TRUE
USE_OMP    = FALSE

USE_REACT  = TRUE

# Include support for participating or tracer particles
USE_PARTICLES = TRUE

SPRAY_COMPONENTS = 8

DEFINES += -DSPRAY_COMPONENTS=$(SPRAY_COMPONENTS)

# This sets the EOS directory in $(PELE_PHYSICS_HOME)/Eos
Eos_dir := Fuego
#Eos_dir     := GammaLaw

# This sets the network directory in $(PELE_PHYSICS_HOME)/Reactions
Reactions_dir := Fuego
#Reactions_dir := Null

# This sets the model directory in $(PELE_PHYSICS_HOME)/Support/Fuego/Mechanism/Models
Chemistry_Model := dodecane_lu

# This sets the transport directory in $(PELE_PHYSICS_HOME)/Transport
Transport_dir := EGLib
#Transport_dir := Constant

Bpack   := ./Make.package 
Blocs   := .

include ${PELEC_HOME}/Exec/Make.PeleC

ifeq ($(USE_PARTICLES), TRUE)
DEFINES += -DPARTICLES
endif

ifeq ($(USE_MPI), TRUE)
 override XTRALIBS += -lifcore
endif
ifeq ($(DEBUG), TRUE)
DEFINES += -DDEBUG
endif
