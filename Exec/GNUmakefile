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

AMREX_HOME ?= /path/to/amrex

EBASE = mini-Castro

DEFINES += -DCRSEGRNDOMP

# This application only supports 3D.
ifneq ($(DIM),3)
  $(error mini-Castro only supports DIM == 3)
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
