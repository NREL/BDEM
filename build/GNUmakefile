AMREX_HOME ?= ../../amrex
BDEM_HOME = ../

EBASE = bdem

DEBUG	= FALSE

DIM	= 3

COMP    = gnu

USE_MPI   = FALSE
USE_OMP   = FALSE
USE_CUDA  = FALSE
USE_EB    = TRUE

BL_NO_FORT = TRUE

TINY_PROFILE = FALSE
USE_PARTICLES = TRUE

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include $(AMREX_HOME)/Src/Base/Make.package
include $(AMREX_HOME)/Src/Particle/Make.package
include $(AMREX_HOME)/Src/EB/Make.package
include $(AMREX_HOME)/Src/Boundary/Make.package
include $(AMREX_HOME)/Src/AmrCore/Make.package

include $(BDEM_HOME)/Source/Make.package
include $(BDEM_HOME)/Source/stl_tools/Make.package
include $(BDEM_HOME)/Source/chemistry/Make.package

INCLUDE_LOCATIONS += $(BDEM_HOME)/Source
INCLUDE_LOCATIONS += $(BDEM_HOME)/Source/stl_tools
INCLUDE_LOCATIONS += $(BDEM_HOME)/Source/chemistry
VPATH_LOCATIONS   += $(BDEM_HOME)/Source
VPATH_LOCATIONS   += $(BDEM_HOME)/Source/stl_tools
VPATH_LOCATIONS   += $(BDEM_HOME)/Source/chemistry

include $(AMREX_HOME)/Tools/GNUMake/Make.rules
