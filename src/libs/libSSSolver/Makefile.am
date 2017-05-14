noinst_LTLIBRARIES = libSSSolver.la
pkginclude_HEADERS = SSSolver.h

libSSSolver_la_SOURCES = 	\
 SSSolver.cc			\
 GetPPI.cc		\
 GetPhiE.cc		\
 StaticSubstrate.cc	\
 SSSolver.h		\
 StaticSubstrate.h	

AM_CPPFLAGS = -I$(top_srcdir)/src/libs/libhrutil    \
              -I$(top_srcdir)/src/libs/libhmat      \
              -I$(top_srcdir)/src/libs/libMatProp   \
              -I$(top_srcdir)/src/libs/libSpherical \
              -I$(top_srcdir)/src/libs/libSGJC      \
              -I$(top_srcdir)/src/libs/libTriInt    \
              -I$(top_srcdir)/src/libs/libMDInterp  \
              -I$(top_srcdir)/src/libs/libscuff     \
              -I$(top_srcdir)/src/libs/libIncField  \
              -I$(top_builddir) # for config.h
