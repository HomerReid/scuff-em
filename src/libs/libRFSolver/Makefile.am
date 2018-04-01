noinst_LTLIBRARIES = libRFSolver.la
pkginclude_HEADERS = RFSolver.h

libRFSolver_la_SOURCES = \
 RWGPorts.cc      \
 MOIIntegrals.cc  \
 ZSConvert.cc     \
 OutputModules.cc


AM_CPPFLAGS = -I$(top_srcdir)/src/libs/libhrutil    \
              -I$(top_srcdir)/src/libs/libhmat      \
              -I$(top_srcdir)/src/libs/libMatProp   \
              -I$(top_srcdir)/src/libs/libSpherical \
              -I$(top_srcdir)/src/libs/libSGJC      \
              -I$(top_srcdir)/src/libs/libTriInt    \
              -I$(top_srcdir)/src/libs/libMDInterp  \
              -I$(top_srcdir)/src/libs/libSubstrate \
              -I$(top_srcdir)/src/libs/libscuff     \
              -I$(top_srcdir)/src/libs/libIncField  \
              -I$(top_builddir) # for config.h
