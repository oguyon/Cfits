#SUFFIXES = .c .h .y .l


bin_PROGRAMS = cfits
cfits_SOURCES = Cfits.c Cfits.h calc_bison.y calc_flex.l calc.h
cfits_CPPFLAGS = -I@abs_top_srcdir@/src -DCONFIGDIR=\"@abs_top_srcdir@/config\" -DSOURCEDIR=\"@abs_top_srcdir@/src\" -DABSSRCTOPDIR=\"@abs_top_srcdir@\" -fopenmp

ACLOCAL_AMFLAGS = -I m4

## For Bison and FLEX
AM_YFLAGS = -d
AM_LFLAGS = -i 

## The contents of this variable are passed to every compilation that invokes the C preprocessor IF target_CPPFLAGS does not exit; it is a list of arguments to the preprocessor. For instance, -I and -D options should be listed here.
## AM_CPPFLAGS = -DCONFIGDIR=\"@abs_top_srcdir@/config\"

noinst_HEADERS = calc.h

BUILT_SOURCES = calc_bison.h


SUBDIRS = 
cfits_LDADD = 


## Modules are built as libraries of convenience
## Linking order matters
## Core modules last, higher level modules at the top

SUBDIRS += AOsystSim
cfits_LDADD += AOsystSim/libaosystsim.la
cfits_CPPFLAGS += -I@abs_top_srcdir@/src/AOsystSim


SUBDIRS += image_format
cfits_LDADD += image_format/libimageformat.la
cfits_CPPFLAGS += -I@abs_top_srcdir@/src/image_format


SUBDIRS += PIAACMCsimul
cfits_LDADD += PIAACMCsimul/libpiaacmcsimul.la
cfits_CPPFLAGS += -I@abs_top_srcdir@/src/PIAACMCsimul

SUBDIRS += coronagraphs
cfits_LDADD += coronagraphs/libcoronagraphs.la
cfits_CPPFLAGS += -I@abs_top_srcdir@/src/coronagraphs

SUBDIRS += ZernikePolyn
cfits_LDADD += ZernikePolyn/libzernikepolyn.la
cfits_CPPFLAGS += -I@abs_top_srcdir@/src/ZernikePolyn

SUBDIRS += WFpropagate
cfits_LDADD += WFpropagate/libwfpropagate.la
cfits_CPPFLAGS += -I@abs_top_srcdir@/src/WFpropagate

SUBDIRS += image_gen
cfits_LDADD += image_gen/libimagegen.la
cfits_CPPFLAGS += -I@abs_top_srcdir@/src/image_gen

SUBDIRS += statistic
cfits_LDADD += statistic/libstatistic.la
cfits_CPPFLAGS += -I@abs_top_srcdir@/src/statistic

SUBDIRS += fft
cfits_LDADD += fft/libfft.la
cfits_CPPFLAGS += -I@abs_top_srcdir@/src/fft

SUBDIRS += info
cfits_LDADD += info/libinfo.la
cfits_CPPFLAGS += -I@abs_top_srcdir@/src/info








## 
## CORE modules, should always be included
##



SUBDIRS += COREMOD_arith
cfits_LDADD += COREMOD_arith/libcoremodarith.la
cfits_CPPFLAGS += -I@abs_top_srcdir@/src/COREMOD_arith

SUBDIRS += COREMOD_iofits
cfits_LDADD += COREMOD_iofits/libcoremodiofits.la
cfits_CPPFLAGS += -I@abs_top_srcdir@/src/COREMOD_iofits

SUBDIRS += COREMOD_memory
cfits_LDADD += COREMOD_memory/libcoremodmemory.la
cfits_CPPFLAGS += -I@abs_top_srcdir@/src/COREMOD_memory

SUBDIRS += COREMOD_tools
cfits_LDADD += COREMOD_tools/libcoremodtools.la
cfits_CPPFLAGS += -I@abs_top_srcdir@/src/COREMOD_tools

SUBDIRS += 00CORE
cfits_LDADD += 00CORE/lib00core.la
cfits_CPPFLAGS += -I@abs_top_srcdir@/src/00CORE

