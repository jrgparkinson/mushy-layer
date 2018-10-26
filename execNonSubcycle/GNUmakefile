# -*- Mode: Makefile -*- 

### This makefile produces an executable for each 
### name in the `ebase' variable

ebase := mushyLayer

#set compiler and CHOMBO_HOME

MACHINE = $(shell uname)
UNAMEM = $(shell uname -m)

DIM=2

TRAPFPE=FALSE
ifeq ($(TRAPFPE),TRUE)	
# petermc, 12 Oct 2007:  never set XTRACPPFLAGS in makefile
# XTRACPPFLAGS += -DTRAP_FPE
endif

# petermc, 12 Oct 2007:  never set XTRACPPFLAGS in makefile
# XTRACPPFLAGS += #-DPAPI 
PAPILIBS = #/usr/local/lib/libpapi.a
PAPILIBS = #/usr/common/usg/papi/2.1/lib/libpapi.a -L/usr/lpp/pmtoolkit/lib -lpmapi

#  have a default for machines not defined here.
# make this a relative path so it will find the Chombo 
# on same level as TDGL
CHOMBO_HOME = ../../../Chombo/lib

include $(CHOMBO_HOME)/mk/Make.defs
include $(CHOMBO_HOME)/mk/Make.defs.config
XTRALIBFLAGS += $(PAPILIBS)

#cases := $(shell echo case*)

##
## names of Chombo libraries needed by this program, in order of search.
##
LibNames :=  AMRElliptic AMRTimeDependent AMRTools BoxTools #OldAMRElliptic  

base_dir = .

TDGL_BASE_DIR = ../

EXEC_DIR = .

ADDEDSRC_DIR = ../srcNonSubcycle
src_dirs = $(ADDEDSRC_DIR)

SRC_DIR = ../src
src_dirs += $(SRC_DIR)

UTIL_DIR = ../util
src_dirs += $(UTIL_DIR)

BCUTIL_DIR = ../BCutil
src_dirs += $(BCUTIL_DIR)



include $(CHOMBO_HOME)/mk/Make.example

ifeq ($(DEBUG), TRUE)
# petermc, 12 Oct 2007:  never set XTRACPPFLAGS in makefile
#XTRACPPFLAGS += -DDEBUG
endif
