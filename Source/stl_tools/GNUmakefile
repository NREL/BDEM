AMREX_HOME ?= ../../amrex

EBASE = stltest

DEBUG	= FALSE

DIM	= 3

COMP    = gcc

USE_MPI   = TRUE
USE_OMP   = FALSE
USE_CUDA  = FALSE
USE_EB    = TRUE

TINY_PROFILE = FALSE

include ./Make.package
CEXE_sources += driver.cpp

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include $(AMREX_HOME)/Src/Base/Make.package
include $(AMREX_HOME)/Src/EB/Make.package
include $(AMREX_HOME)/Src/Boundary/Make.package
include $(AMREX_HOME)/Src/AmrCore/Make.package

include $(AMREX_HOME)/Tools/GNUMake/Make.rules
