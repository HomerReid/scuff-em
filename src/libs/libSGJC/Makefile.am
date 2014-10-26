noinst_LTLIBRARIES = libSGJC.la
pkginclude_HEADERS = libSGJC.h
libSGJC_la_SOURCES = \
 pcubature.c \
 hcubature.c \
 libSGJC.c   \
 clencurt.h  \
 converged.h \
 vwrapper.h

EXTRA_DIST = COPYING README

noinst_PROGRAMS = tlibSGJC tlibSGJC2 # tlibVVQAG

tlibSGJC_SOURCES = tlibSGJC.cc
tlibSGJC_LDADD = libSGJC.la
tlibSGJC2_SOURCES = tlibSGJC2.cc
tlibSGJC2_LDADD = libSGJC.la
# tlibVVQAG_SOURCES = tlibVVQAG.cc
# tlibVVQAG_LDADD = libSGJC.la
