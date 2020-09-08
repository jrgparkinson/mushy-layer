# -*- Mode: Makefile -*- 

## Define the variables needed by Make.example

# the base name(s) of the application(s) in this directory
ebase = mushyLayer

# a bunch of common defs for the various machines we might be running on
MACHINE := $(shell uname)
UNAMEM := $(shell uname -m)
UNAMEN := $(shell uname -n)
WHO := $(shell whoami)

# the location of the Chombo "lib" directory
ifeq ($(NERSC_HOST),cori) 
CHOMBO_HOME = ../../Chombo/lib
else
CHOMBO_HOME = ../../chombofork/lib
endif

ifeq ($(UNAMEN),martin)
# Compile with chombo trunk to check everything will work for other users
CHOMBO_HOME = ../../Chombo/lib
endif


include ../mk/Make.defs.MushyLayer


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

# CUSTOM_CHOMBO_LIB=true

# application-specific variables

# application-specific targets
