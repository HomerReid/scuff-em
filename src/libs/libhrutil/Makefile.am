noinst_LTLIBRARIES = libhrutil.la
pkginclude_HEADERS = libhrutil.h
libhrutil_la_SOURCES = \
 libhrutil.cc         \
 ProcessArguments.cc  \
 ProcessOptions.cc    \
 Vector.cc

noinst_PROGRAMS = tProcessArguments tProcessOptions
tProcessArguments_SOURCES = tProcessArguments.cc
tProcessArguments_LDADD = libhrutil.la
tProcessOptions_SOURCES = tProcessOptions.cc
tProcessOptions_LDADD = libhrutil.la

AM_CPPFLAGS = -I$(top_builddir) # for config.h 
