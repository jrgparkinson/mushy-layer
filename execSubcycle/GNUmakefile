# -*- Mode: Makefile -*- 

# trace the chain of included makefiles
#makefiles += releasedExamples_AMRGodunov_execAdvectDiffuse

## Define the variables needed by Make.example

# the base name(s) of the application(s) in this directory
ebase = mushyLayer

# the location of the Chombo "lib" directory
CHOMBO_HOME = ../../chombofork/lib

include $(CHOMBO_HOME)/mk/Make.defs
include $(CHOMBO_HOME)/mk/Make.defs.config

# names of Chombo libraries needed by this program, in order of search.
LibNames = AMRElliptic AMRTimeDependent AMRTools BoxTools

# the locations of the source code directories
EXEC_DIR = .
base_dir = .
src_dirs = ../srcSubcycle ../util ../src  ../BCutil 
# ../muparser-2.2.5/include


# input file for 'run' target
INPUT = inputs

# shared code for building example programs
include $(CHOMBO_HOME)/mk/Make.example

# application-specific variables

# application-specific targets
