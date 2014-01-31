bin_PROGRAMS = scuff-neq

scuff_neq_SOURCES = 		\
 GetFlux.cc             	\
 CreateSNEQData.cc       	\
 SIPFT.cc              		\
 Quadrature.cc         		\
 scuff-neq.cc         		\
 scuff-neq.h

scuff_neq_LDADD = $(top_builddir)/src/libs/libscuff/libscuff.la

AM_CPPFLAGS = -I$(top_srcdir)/src/libs/libscuff      \
              -I$(top_srcdir)/src/libs/libIncField   \
              -I$(top_srcdir)/src/libs/libMatProp    \
              -I$(top_srcdir)/src/libs/libMDInterp   \
              -I$(top_srcdir)/src/libs/libhmat       \
              -I$(top_srcdir)/src/libs/libSGJC       \
              -I$(top_srcdir)/src/libs/libTriInt     \
              -I$(top_srcdir)/src/libs/libhrutil
