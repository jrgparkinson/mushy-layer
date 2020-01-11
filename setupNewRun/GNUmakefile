# -*- Mode: Makefile -*- 

### This makefile produces an executable for each 
### name in the `ebase' variable
ebase = setupnewrun

## location of 'lib' directory
#CHOMBO_HOME = ../../chombofork/lib


include ../mk/Make.defs.MushyLayer


##
## names of Chombo libraries needed by this program, in order of search.
##
LibNames = AMRElliptic AMRTimeDependent AMRTools BoxTools

#LibNames = AMRTools BoxTools

EXEC_DIR = .
base_dir = .
src_dirs = ../srcSubcycle/ ../util ../src ../BCutil

INPUT = inputs

include $(CHOMBO_HOME)/mk/Make.test


override VERBOSE = # #program doesnt handle -q option



