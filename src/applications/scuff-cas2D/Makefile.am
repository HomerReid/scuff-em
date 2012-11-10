SUBDIRS = libTDRT libAmosBessel

bin_PROGRAMS = scuff-cas2D

scuff_cas2D_SOURCES = 		\
 scuff-cas2D.h  		\
 scuff-cas2D.cc			\
 CreateC2DWorkspace.cc		\
 CubatureMethods.cc		\
 ImageObjects.cc		\
 XQIntegrand.cc

scuff_cas2D_LDADD = \
 $(builddir)/libTDRT/libTDRT.la				\
 $(builddir)/libAmosBessel/libAmosBessel.la		\
 $(top_builddir)/src/libs/libSGJC/libSGJC.la		\
 $(top_builddir)/src/libs/libMDInterp/libMDInterp.la	\
 $(top_builddir)/src/libs/libMatProp/libMatProp.la	\
 $(top_builddir)/src/libs/libhmat/libhmat.la		\
 $(top_builddir)/src/libs/libhrutil/libhrutil.la


AM_CPPFLAGS = \
  -I$(top_srcdir)/src/libs/libhrutil				\
  -I$(top_srcdir)/src/applications/scuff-cas2D/libTDRT		\
  -I$(top_srcdir)/src/applications/scuff-cas2D/libAmosBessel 	\
  -I$(top_srcdir)/src/libs/libSGJC       \
  -I$(top_srcdir)/src/libs/libMDInterp   \
  -I$(top_srcdir)/src/libs/libMatProp    \
  -I$(top_srcdir)/src/libs/libhmat
