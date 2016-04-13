bin_PROGRAMS = scuff-neq scuff-integrate

scuff_neq_SOURCES = 		\
 scuff-neq.cc         		\
 scuff-neq.h			\
 GetFlux.cc             	\
 CreateSNEQData.cc       	\
 Quadrature.cc

scuff_neq_LDADD = $(top_builddir)/src/libs/libscuff/libscuff.la

scuff_integrate_SOURCES = scuff-integrate.cc

scuff_integrate_LDADD = $(top_builddir)/src/libs/libscuff/libscuff.la

AM_CPPFLAGS = -I$(top_srcdir)/src/libs/libscuff      \
              -I$(top_srcdir)/src/libs/libIncField   \
              -I$(top_srcdir)/src/libs/libMatProp    \
              -I$(top_srcdir)/src/libs/libMDInterp   \
              -I$(top_srcdir)/src/libs/libhmat       \
              -I$(top_srcdir)/src/libs/libSGJC       \
              -I$(top_srcdir)/src/libs/libTriInt     \
              -I$(top_srcdir)/src/libs/libhrutil
