bin_PROGRAMS = scuff-caspol

scuff_caspol_SOURCES = 		\
 scuff-caspol.h			\
 scuff-caspol.cc          	\
 GetCPIntegrand.cc 		\
 PolModel.cc
scuff_caspol_LDADD = $(top_builddir)/src/libs/libscuff/libscuff.la

EXTRA_DIST = PolData.dat

AM_CPPFLAGS = -I$(top_srcdir)/src/libs/libscuff        \
              -I$(top_srcdir)/src/libs/libIncField   \
              -I$(top_srcdir)/src/libs/libMatProp    \
              -I$(top_srcdir)/src/libs/libMDInterp    \
              -I$(top_srcdir)/src/libs/libhmat       \
              -I$(top_srcdir)/src/libs/libSGJC       \
              -I$(top_srcdir)/src/libs/libTriInt     \
              -I$(top_srcdir)/src/libs/libMDInterp   \
              -I$(top_srcdir)/src/libs/libhrutil
