bin_PROGRAMS = cfits
cfits_SOURCES = Cfits.c Cfits.h
cfits_CPPFLAGS = -I@abs_top_srcdir@/src

AM_CPPFLAGS = -DCONFIGDIR=\"tttest\"


SUBDIRS = 00CORE
#cfits_SOURCES += 00CORE/00CORE.c 00CORE/00CORE.h
cfits_LDADD = 00CORE/lib00core.la
cfits_CPPFLAGS += -I@abs_top_srcdir@/src/00CORE

SUBDIRS += COREMOD_arith
#cfits_SOURCES += COREMOD_arith/COREMOD_arith.c COREMOD_arith/COREMOD_arith.h
cfits_LDADD += COREMOD_arith/libcoremodarith.la
cfits_CPPFLAGS += -I@abs_top_srcdir@/src/COREMOD_arith

SUBDIRS += COREMOD_fft
#cfits_SOURCES += COREMOD_fft/COREMOD_fft.c COREMOD_fft/COREMOD_fft.h
cfits_LDADD += COREMOD_fft/libcoremodfft.la
cfits_CPPFLAGS += -I@abs_top_srcdir@/src/COREMOD_fft

SUBDIRS += COREMOD_info
#cfits_SOURCES += COREMOD_info/COREMOD_info.c COREMOD_info/COREMOD_info.h
cfits_LDADD += COREMOD_info/libcoremodinfo.la
cfits_CPPFLAGS += -I@abs_top_srcdir@/src/COREMOD_info

SUBDIRS += COREMOD_iofits
#cfits_SOURCES += COREMOD_iofits/COREMOD_iofits.c COREMOD_iofits/COREMOD_iofits.h
cfits_LDADD += COREMOD_iofits/libcoremodiofits.la
cfits_CPPFLAGS += -I@abs_top_srcdir@/src/COREMOD_iofits

SUBDIRS += COREMOD_memory
#cfits_SOURCES += COREMOD_memory/COREMOD_memory.c COREMOD_memory/COREMOD_memory.h
cfits_LDADD += COREMOD_memory/libcoremodmemory.la
cfits_CPPFLAGS += -I@abs_top_srcdir@/src/COREMOD_memory

SUBDIRS += COREMOD_tools
#cfits_SOURCES += COREMOD_tools/COREMOD_tools.c COREMOD_tools/COREMOD_tools.h
cfits_LDADD += COREMOD_tools/libcoremodtools.la
cfits_CPPFLAGS += -I@abs_top_srcdir@/src/COREMOD_tools

