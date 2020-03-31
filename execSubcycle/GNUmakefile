# -*- Mode: Makefile -*- 

## Define the variables needed by Make.example

# the base name(s) of the application(s) in this directory
ebase = mushyLayer

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
