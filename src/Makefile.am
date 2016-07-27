#SUFFIXES = .c .h .y .l

bin_PROGRAMS = cfitsTK
cfitsTK_SOURCES = CLIcore.c CLIcore.h initmodules.c initmodules.h calc_bison.y calc_flex.l calc.h
cfitsTK_CPPFLAGS = -I@abs_top_srcdir@/src -DCONFIGDIR=\"@abs_top_srcdir@/config\" -DSOURCEDIR=\"@abs_top_srcdir@/src\" -DABSSRCTOPDIR=\"@abs_top_srcdir@\" -fopenmp

ACLOCAL_AMFLAGS = -I m4

## For Bison and FLEX
AM_YFLAGS = -d
AM_LFLAGS = -i 

## The contents of this variable are passed to every compilation that invokes the C preprocessor IF target_CPPFLAGS does not exit; it is a list of arguments to the preprocessor. For instance, -I and -D options should be listed here.
## AM_CPPFLAGS = -DCONFIGDIR=\"@abs_top_srcdir@/config\"

noinst_HEADERS = calc.h

BUILT_SOURCES = calc_bison.h


SUBDIRS = 
cfitsTK_LDADD = 


## Modules are built as libraries of convenience
## Linking order matters
## Core modules last, higher level modules at the top


SUBDIRS += linARfilterPred 
cfitsTK_LDADD += linARfilterPred/liblinarfilterpred.la
cfitsTK_CPPFLAGS += -I@abs_top_srcdir@/src/linARfilterPred

SUBDIRS += TransitLC 
cfitsTK_LDADD += TransitLC/libtransitlc.la
cfitsTK_CPPFLAGS += -I@abs_top_srcdir@/src/TransitLC

SUBDIRS += SCExAO_control
cfitsTK_LDADD += SCExAO_control/libscexaocontrol.la
cfitsTK_CPPFLAGS += -I@abs_top_srcdir@/src/SCExAO_control

SUBDIRS += cudacomp
cfitsTK_LDADD += cudacomp/libcudacomp.la
cfitsTK_CPPFLAGS += -I@abs_top_srcdir@/src/cudacomp

SUBDIRS += AtmosphericTurbulence
cfitsTK_LDADD += AtmosphericTurbulence/libatmosphericturbulence.la
cfitsTK_CPPFLAGS += -I@abs_top_srcdir@/src/AtmosphericTurbulence

SUBDIRS += AtmosphereModel
cfitsTK_LDADD += AtmosphereModel/libatmospheremodel.la
cfitsTK_CPPFLAGS += -I@abs_top_srcdir@/src/AtmosphereModel


SUBDIRS += psf
cfitsTK_LDADD += psf/libpsf.la
cfitsTK_CPPFLAGS += -I@abs_top_srcdir@/src/psf

SUBDIRS += AOloopControl
cfitsTK_LDADD += AOloopControl/libaoloopcontrol.la
cfitsTK_CPPFLAGS += -I@abs_top_srcdir@/src/AOloopControl

SUBDIRS += AOsystSim
cfitsTK_LDADD += AOsystSim/libaosystsim.la
cfitsTK_CPPFLAGS += -I@abs_top_srcdir@/src/AOsystSim

SUBDIRS += AOloopControl_DM
cfitsTK_LDADD += AOloopControl_DM/libaoloopcontroldm.la
cfitsTK_CPPFLAGS += -I@abs_top_srcdir@/src/AOloopControl_DM

SUBDIRS += img_reduce
cfitsTK_LDADD += img_reduce/libimgreduce.la
cfitsTK_CPPFLAGS += -I@abs_top_srcdir@/src/img_reduce

SUBDIRS += image_format
cfitsTK_LDADD += image_format/libimageformat.la
cfitsTK_CPPFLAGS += -I@abs_top_srcdir@/src/image_format

SUBDIRS += PIAACMCsimul
cfitsTK_LDADD += PIAACMCsimul/libpiaacmcsimul.la
cfitsTK_CPPFLAGS += -I@abs_top_srcdir@/src/PIAACMCsimul

SUBDIRS += OptSystProp
cfitsTK_LDADD += OptSystProp/liboptsystprop.la
cfitsTK_CPPFLAGS += -I@abs_top_srcdir@/src/OptSystProp

SUBDIRS += OpticsMaterials
cfitsTK_LDADD += OpticsMaterials/libopticsmaterials.la
cfitsTK_CPPFLAGS += -I@abs_top_srcdir@/src/OpticsMaterials

SUBDIRS += coronagraphs
cfitsTK_LDADD += coronagraphs/libcoronagraphs.la
cfitsTK_CPPFLAGS += -I@abs_top_srcdir@/src/coronagraphs

SUBDIRS += ZernikePolyn
cfitsTK_LDADD += ZernikePolyn/libzernikepolyn.la
cfitsTK_CPPFLAGS += -I@abs_top_srcdir@/src/ZernikePolyn

SUBDIRS += WFpropagate
cfitsTK_LDADD += WFpropagate/libwfpropagate.la
cfitsTK_CPPFLAGS += -I@abs_top_srcdir@/src/WFpropagate

SUBDIRS += image_basic
cfitsTK_LDADD += image_basic/libimagebasic.la
cfitsTK_CPPFLAGS += -I@abs_top_srcdir@/src/image_basic

SUBDIRS += image_filter
cfitsTK_LDADD += image_filter/libimagefilter.la
cfitsTK_CPPFLAGS += -I@abs_top_srcdir@/src/image_filter

SUBDIRS += kdtree
cfitsTK_LDADD += kdtree/libkdtree.la
cfitsTK_CPPFLAGS += -I@abs_top_srcdir@/src/kdtree

SUBDIRS += image_gen
cfitsTK_LDADD += image_gen/libimagegen.la
cfitsTK_CPPFLAGS += -I@abs_top_srcdir@/src/image_gen

SUBDIRS += linopt_imtools
cfitsTK_LDADD += linopt_imtools/liblinoptimtools.la
cfitsTK_CPPFLAGS += -I@abs_top_srcdir@/src/linopt_imtools

SUBDIRS += statistic
cfitsTK_LDADD += statistic/libstatistic.la
cfitsTK_CPPFLAGS += -I@abs_top_srcdir@/src/statistic

SUBDIRS += fft
cfitsTK_LDADD += fft/libfft.la
cfitsTK_CPPFLAGS += -I@abs_top_srcdir@/src/fft

SUBDIRS += info
cfitsTK_LDADD += info/libinfo.la
cfitsTK_CPPFLAGS += -I@abs_top_srcdir@/src/info








## 
## CORE modules, should always be included
##



SUBDIRS += COREMOD_arith
cfitsTK_LDADD += COREMOD_arith/libcoremodarith.la
cfitsTK_CPPFLAGS += -I@abs_top_srcdir@/src/COREMOD_arith

SUBDIRS += COREMOD_iofits
cfitsTK_LDADD += COREMOD_iofits/libcoremodiofits.la
cfitsTK_CPPFLAGS += -I@abs_top_srcdir@/src/COREMOD_iofits

SUBDIRS += COREMOD_memory
cfitsTK_LDADD += COREMOD_memory/libcoremodmemory.la
cfitsTK_CPPFLAGS += -I@abs_top_srcdir@/src/COREMOD_memory

SUBDIRS += COREMOD_tools
cfitsTK_LDADD += COREMOD_tools/libcoremodtools.la
cfitsTK_CPPFLAGS += -I@abs_top_srcdir@/src/COREMOD_tools

SUBDIRS += 00CORE
cfitsTK_LDADD += 00CORE/lib00core.la
cfitsTK_CPPFLAGS += -I@abs_top_srcdir@/src/00CORE

