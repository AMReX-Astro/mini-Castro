PRECISION  = DOUBLE
PROFILE    = FALSE

DEBUG      = FALSE

DIM        = 3

COMP	   = PGI

USE_MPI    = FALSE
USE_OMP    = FALSE
USE_CUDA   = FALSE

TINY_PROFILE = FALSE

USE_GPU_PRAGMA = TRUE

# define the location of the CASTRO top directory
CASTRO_HOME  := ../
TOP := $(CASTRO_HOME)

AMREX_HOME ?= /path/to/amrex

EBASE = Castro

# This application only supports 3D.
ifneq ($(DIM),3)
  $(error mini-Castro only supports BL_SPACEDIM == 3)
endif

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

# Castro uses a coarse grained OMP approach
DEFINES += -DCRSEGRNDOMP

ifeq ($(USE_ACC), TRUE)
   DEFINES += -DACC
endif

all: $(executable)
	@echo SUCCESS

# AMReX directories -- note the Make.package for these adds these
# directories into VPATH_LOCATIONS and INCLUDE_LOCATIONS for us, so we
# don't need to do it here
Pdirs 	:= Base AmrCore Amr Boundary

Bpack	+= $(foreach dir, $(Pdirs), $(AMREX_HOME)/Src/$(dir)/Make.package)

# Castro directories
Bdirs 	:= Source

Bpack	+= $(foreach dir, $(Bdirs), $(TOP)/$(dir)/Make.package)
Blocs	+= $(foreach dir, $(Bdirs), $(TOP)/$(dir))

include $(Bpack)

INCLUDE_LOCATIONS += $(Blocs)
VPATH_LOCATIONS   += $(Blocs)

include $(AMREX_HOME)/Tools/GNUMake/Make.rules
