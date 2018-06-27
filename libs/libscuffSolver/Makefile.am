noinst_LTLIBRARIES = libscuffSolver.la
pkginclude_HEADERS = scuffSolver.h

libscuffSolver_la_SOURCES = 	\
 scuffSolver.h            	\
 MOIIntegrals.cc  		\
 OutputModules.cc		\
 scuffSolver.cc      		\
 RWGPorts.cc      		

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

AM_CXXFLAGS = -std=c++11

SUBDIRS = tests
