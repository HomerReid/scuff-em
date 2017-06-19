bin_PROGRAMS = scuff-static

scuff_static_SOURCES = 		\
 scuff-static.h			\
 BMAccelerator.cc		\
 OutputModules.cc		\
 scuff-static.cc

scuff_static_LDADD = $(top_builddir)/src/libs/libscuff/libscuff.la

AM_CPPFLAGS = -I$(top_srcdir)/src/libs/libscuff      \
              -I$(top_srcdir)/src/libs/libIncField   \
              -I$(top_srcdir)/src/libs/libMatProp    \
              -I$(top_srcdir)/src/libs/libMatProp/cmatheval \
              -I$(top_srcdir)/src/libs/libMDInterp   \
              -I$(top_srcdir)/src/libs/libhmat       \
              -I$(top_srcdir)/src/libs/libSGJC       \
              -I$(top_srcdir)/src/libs/libSpherical  \
              -I$(top_srcdir)/src/libs/libTriInt     \
              -I$(top_srcdir)/src/libs/libSSSolver   \
              -I$(top_srcdir)/src/libs/libhrutil
