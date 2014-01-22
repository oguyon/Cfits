bin_PROGRAMS = cfits
cfits_SOURCES = Cfits.c Cfits.h
cfits_CPPFLAGS = -I@abs_top_srcdir@/src

AM_CPPFLAGS = -DCONFIGDIR=\"@abs_top_srcdir@/config\"


SUBDIRS = 00CORE
cfits_LDADD = 00CORE/lib00core.la
cfits_CPPFLAGS += -I@abs_top_srcdir@/src/00CORE

SUBDIRS += COREMOD_fft
cfits_LDADD += COREMOD_fft/libcoremodfft.la
cfits_CPPFLAGS += -I@abs_top_srcdir@/src/COREMOD_fft

SUBDIRS += COREMOD_arith
cfits_LDADD += COREMOD_arith/libcoremodarith.la
cfits_CPPFLAGS += -I@abs_top_srcdir@/src/COREMOD_arith

SUBDIRS += COREMOD_info
cfits_LDADD += COREMOD_info/libcoremodinfo.la
cfits_CPPFLAGS += -I@abs_top_srcdir@/src/COREMOD_info

SUBDIRS += COREMOD_iofits
cfits_LDADD += COREMOD_iofits/libcoremodiofits.la
cfits_CPPFLAGS += -I@abs_top_srcdir@/src/COREMOD_iofits

SUBDIRS += COREMOD_memory
cfits_LDADD += COREMOD_memory/libcoremodmemory.la
cfits_CPPFLAGS += -I@abs_top_srcdir@/src/COREMOD_memory

SUBDIRS += COREMOD_tools
cfits_LDADD += COREMOD_tools/libcoremodtools.la
cfits_CPPFLAGS += -I@abs_top_srcdir@/src/COREMOD_tools

