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

TINY_PROFILE = FALSE

USE_GPU_PRAGMA = TRUE

EBASE = mini-Castro

DEFINES += -DCRSEGRNDOMP

AMREX_HOME := ../amrex

# This application only supports 3D.
ifneq ($(DIM),3)
  $(error mini-Castro only supports DIM == 3)
endif

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
