PRECISION  = DOUBLE
PROFILE    = FALSE

DEBUG      = FALSE

DIM        = 3

COMP	   = PGI

USE_MPI    = FALSE
USE_OMP    = FALSE
USE_CUDA   = FALSE

TINY_PROFILE = FALSE

# define the location of the CASTRO top directory
CASTRO_HOME  := ../

AMREX_HOME ?= /path/to/amrex
CASTRO_HOME ?= /path/to/Castro

TOP := $(CASTRO_HOME)

EBASE = Castro

# This application only supports 3D.
ifneq ($(DIM),3)
  $(error StarLord only supports BL_SPACEDIM == 3)
endif

# EXTERN_CORE is simply a list of the directories we use that
# will be searched for runtime parameters
EXTERN_CORE ?= $(TOP)/Source

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

# Castro uses a coarse grained OMP approach
DEFINES += -DCRSEGRNDOMP

# Limit the maximum number of registers available
ifeq ($(lowercase_comp),pgi)
  ifeq ($(USE_CUDA),TRUE)
    FFLAGS += -Mcuda=maxregcount:128
    F90FLAGS += -Mcuda=maxregcount:128
  endif
endif

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

ifeq ($(MEGA_FORTRAN), TRUE)
   F90EXE_sources += mega_f.F90
else
   F90EXE_sources += $(ca_F90EXE_sources)
endif

mega_f.F90: $(ca_F90EXE_sources)
	$(AMREX_HOME)/Tools/F_scripts/f90cat.py $^

include $(AMREX_HOME)/Src/F_BaseLib/FParallelMG.mak

# runtime parameter support for extern/ routines


EXTERN_TEMPLATE := $(TOP)/Source/extern_probin.template

EXTERN_PARAMETERS := $(shell $(AMREX_HOME)/Tools/F_scripts/findparams.py $(EXTERN_CORE))

# job_info support
CEXE_headers += $(AMREX_HOME)/Tools/C_scripts/AMReX_buildInfo.H
INCLUDE_LOCATIONS +=  $(AMREX_HOME)/Tools/C_scripts

# we make AMReX_buildInfo.cpp as we make the .o file, so we can delete it
# immediately.  this way if the build is interrupted, we are guaranteed
# to remake it

objForExecs += $(objEXETempDir)/AMReX_buildInfo.o

.FORCE:
.PHONE: .FORCE

$(objEXETempDir)/AMReX_buildInfo.o: .FORCE
	$(AMREX_HOME)/Tools/C_scripts/makebuildinfo_C.py \
          --amrex_home "$(AMREX_HOME)" \
          --COMP "$(COMP)" --COMP_VERSION "$(COMP_VERSION)" \
          --CXX_comp_name "$(CXX)" --CXX_flags "$(CXXFLAGS) $(CPPFLAGS) $(includes)" \
          --F_comp_name "$(F90)" --F_flags "$(F90FLAGS)" \
          --link_flags "$(LDFLAGS)" --libraries "$(libraries)" \
          --MODULES "$(MNAMES)" $(EXTRA_BUILD_INFO) \
          --GIT "$(TOP) $(AMREX_HOME)"
	$(SILENT) $(CCACHE) $(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $(CXXEXEFLAGS) AMReX_buildInfo.cpp -o $(objEXETempDir)/AMReX_buildInfo.o
	$(SILENT) $(RM) AMReX_buildInfo.cpp


include $(AMREX_HOME)/Tools/GNUMake/Make.rules

clean::
	$(SILENT) $(RM) AMReX_buildInfo.cpp

#-----------------------------------------------------------------------------
# for debugging.  To see the value of a Makefile variable,
# e.g. Fmlocs, simply do "make print-Fmlocs".  This will
# print out the value.
print-%::
	@echo "$* is $($*)"
	@$(RM) AMReX_buildInfo.cpp
