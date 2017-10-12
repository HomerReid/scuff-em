bin_PROGRAMS = scuff-ldos

scuff_ldos_SOURCES = 		\
 scuff-ldos.h			\
 scuff-ldos.cc         		\
 CreateSLDData.cc 		\
 GetLDOS.cc			\
 AnalyticalDGFs.cc

scuff_ldos_LDADD = $(top_builddir)/src/libs/libscuff/libscuff.la

AM_CPPFLAGS = -I$(top_srcdir)/src/libs/libscuff      \
              -I$(top_srcdir)/src/libs/libIncField   \
              -I$(top_srcdir)/src/libs/libMatProp    \
              -I$(top_srcdir)/src/libs/libMDInterp   \
              -I$(top_srcdir)/src/libs/libhmat       \
              -I$(top_srcdir)/src/libs/libSGJC       \
              -I$(top_srcdir)/src/libs/libSubstrate  \
              -I$(top_srcdir)/src/libs/libSpherical  \
              -I$(top_srcdir)/src/libs/libTriInt     \
              -I$(top_srcdir)/src/libs/libhrutil
