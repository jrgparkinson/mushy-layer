# -*- Mode: Makefile -*- 

## Define the variables needed by Make.example

# the base name(s) of the application(s) in this directory
ebase = mushyLayer

# the location of the Chombo "lib" directory
CHOMBO_HOME = ../../chombofork/lib
# Compile with chombo trunk to check everything will work for other users
#CHOMBO_HOME = ../../Chombo-trunk/lib 

HOST=$(shell hostname)
GYRE=gyre
LOCAL=atmlxlap005

# If variable GYRE is contained in HOST (e.g. HOST=gyre1, gyre2, ...), then add .GYRE to the compiled names
# This is because we can only run code on the gyre servers which has been compiled there
ifneq (,$(findstring $(GYRE),$(HOST)))
    XTRACONFIG = .GYRE
endif

#ifneq (,$(findstring $(LOCAL),$(HOST)))
#    XTRACONFIG = .LOCAL
#endif

#XTRACONFIG= ifeq ($(GCC_MINOR),$(filter $(GCC_MINOR),4 5)) .GYRE


include $(CHOMBO_HOME)/mk/Make.defs
include $(CHOMBO_HOME)/mk/Make.defs.config

# names of Chombo libraries needed by this program, in order of search.
LibNames = AMRElliptic AMRTimeDependent AMRTools BoxTools

# the locations of the source code directories
EXEC_DIR = .
base_dir = .
src_dirs = ../srcSubcycle ../util ../src  ../BCutil 

# input file for 'run' target
INPUT = inputs

# shared code for building example programs
include $(CHOMBO_HOME)/mk/Make.example

CUSTOM_CHOMBO_LIB=true

# application-specific variables

# application-specific targets
