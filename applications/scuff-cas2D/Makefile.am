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
 $(top_builddir)/libs/libSGJC/libSGJC.la		\
 $(top_builddir)/libs/libMDInterp/libMDInterp.la	\
 $(top_builddir)/libs/libMatProp/libMatProp.la	\
 $(top_builddir)/libs/libhmat/libhmat.la		\
 $(top_builddir)/libs/libhrutil/libhrutil.la


AM_CPPFLAGS = \
  -I$(top_srcdir)/libs/libhrutil				\
  -I$(top_srcdir)/applications/scuff-cas2D/libTDRT		\
  -I$(top_srcdir)/applications/scuff-cas2D/libAmosBessel 	\
  -I$(top_srcdir)/libs/libSGJC       \
  -I$(top_srcdir)/libs/libMDInterp   \
  -I$(top_srcdir)/libs/libMatProp    \
  -I$(top_srcdir)/libs/libhmat     \
  -DVERSIONSTRING=\"$(shell git describe always --dirty)\"
