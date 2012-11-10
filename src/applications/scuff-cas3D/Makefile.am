bin_PROGRAMS = scuff-cas3D

scuff_cas3D_SOURCES = 		\
 CasimirIntegrand.cc 		\
 CreateSC3Data.cc       	\
 SumsIntegrals.cc       	\
 scuff-cas3D.cc         	\
 scuff-cas3D.h

scuff_cas3D_LDADD = $(top_builddir)/src/libs/libscuff/libscuff.la

AM_CPPFLAGS = -I$(top_srcdir)/src/libs/libscuff      \
              -I$(top_srcdir)/src/libs/libIncField   \
              -I$(top_srcdir)/src/libs/libMatProp    \
              -I$(top_srcdir)/src/libs/libMDInterp   \
              -I$(top_srcdir)/src/libs/libhmat       \
              -I$(top_srcdir)/src/libs/libSGJC       \
              -I$(top_srcdir)/src/libs/libTriInt     \
              -I$(top_srcdir)/src/libs/libhrutil
