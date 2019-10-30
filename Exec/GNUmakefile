PRECISION  = DOUBLE
PROFILE    = FALSE

DEBUG      = FALSE

DIM        = 3

COMP	   = PGI

# To work around PGI compiler bug
PGI_GOPT   = FALSE

USE_MPI    = FALSE
USE_OMP    = FALSE
USE_CUDA   = FALSE
USE_ACC    = FALSE

CUDA_VERBOSE = FALSE

# We only support OpenACC if CUDA is also defined.
# This is required because AMReX uses CUDA internally
# for its operations, and those would massively slow
# down the code on GPUs if they're not accelerated.

ifeq ($(USE_ACC),TRUE)
  USE_CUDA = TRUE
endif

ifeq ($(USE_ACC),TRUE)
  DEFINES += -DCASTRO_DEVICE=
  DEFINES += -DCASTRO_FORT_DEVICE=
else
  DEFINES += -DCASTRO_DEVICE=AMREX_GPU_DEVICE
  DEFINES += -DCASTRO_FORT_DEVICE=AMREX_CUDA_FORT_DEVICE
endif

TINY_PROFILE = FALSE

EBASE = mini-Castro

DEFINES += -DCRSEGRNDOMP

# This application only supports 3D.
ifneq ($(DIM),3)
  $(error mini-Castro only supports DIM == 3)
endif

# If the user doesn't provide AMReX, use the git submodule version.
AMREX_HOME ?= ../amrex

# Include the AMReX make rules. Throw an error
# if we don't have AMReX.
ifeq ("$(wildcard $(AMREX_HOME)/Tools/GNUMake/Make.defs)","")
  $(error AMReX has not been downloaded. Please run "git submodule update --init" from the top level of the code)
endif
include $(AMREX_HOME)/Tools/GNUMake/Make.defs

all: $(executable)
	@echo SUCCESS

# AMReX directories
Pdirs 	:= Base AmrCore Amr Boundary

Bpack	+= $(foreach dir, $(Pdirs), $(AMREX_HOME)/Src/$(dir)/Make.package)

# mini-Castro directory
Bdirs 	:= Source

Bpack	+= ../Source/Make.package
Blocs	+= ../Source

include $(Bpack)

INCLUDE_LOCATIONS += $(Blocs)
VPATH_LOCATIONS   += $(Blocs)

include $(AMREX_HOME)/Tools/GNUMake/Make.rules
