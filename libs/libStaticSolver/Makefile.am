noinst_LTLIBRARIES = libStaticSolver.la
pkginclude_HEADERS = StaticSolver.h

libStaticSolver_la_SOURCES = 	\
 StaticSolver.cc	\
 Excitations.cc    	\
 GetPPI.cc		\
 GetPhiE.cc		\
 PostProcessing.cc

AM_CPPFLAGS = -I$(top_srcdir)/libs/libhrutil    \
              -I$(top_srcdir)/libs/libhmat      \
              -I$(top_srcdir)/libs/libMatProp   \
              -I$(top_srcdir)/libs/libSpherical \
              -I$(top_srcdir)/libs/libSGJC      \
              -I$(top_srcdir)/libs/libTriInt    \
              -I$(top_srcdir)/libs/libMDInterp  \
              -I$(top_srcdir)/libs/libSubstrate \
              -I$(top_srcdir)/libs/libscuff     \
              -I$(top_srcdir)/libs/libIncField  \
              -I$(top_builddir) # for config.h
