bin_PROGRAMS = scuff-cas3D

scuff_cas3D_SOURCES = 		\
 CasimirIntegrand.cc 		\
 CreateSC3Data.cc       	\
 SumsIntegrals.cc       	\
 scuff-cas3D.cc         	\
 scuff-cas3D.h

scuff_cas3D_LDADD = $(top_builddir)/libs/libscuff/libscuff.la

AM_CPPFLAGS = -I$(top_srcdir)/libs/libscuff      \
              -I$(top_srcdir)/libs/libIncField   \
              -I$(top_srcdir)/libs/libMatProp    \
              -I$(top_srcdir)/libs/libMDInterp   \
              -I$(top_srcdir)/libs/libhmat       \
              -I$(top_srcdir)/libs/libSGJC       \
              -I$(top_srcdir)/libs/libSubstrate  \
              -I$(top_srcdir)/libs/libTriInt     \
              -I$(top_srcdir)/libs/libhrutil     \
              -DVERSIONSTRING=\"$(shell git describe always --dirty)\"
