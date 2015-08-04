bin_PROGRAMS = scuff-ldos

scuff_ldos_SOURCES = 		\
 scuff-ldos.cc         		\
 scuff-ldos.h			\
 CreateSLDData.cc 		\
 GetLDOS.cc			\
 Integration.cc			\
 AnalyticalDGFs.cc		\
 LatticeSum.cc    

scuff_ldos_LDADD = $(top_builddir)/src/libs/libscuff/libscuff.la

AM_CPPFLAGS = -I$(top_srcdir)/src/libs/libscuff      \
              -I$(top_srcdir)/src/libs/libIncField   \
              -I$(top_srcdir)/src/libs/libMatProp    \
              -I$(top_srcdir)/src/libs/libMDInterp   \
              -I$(top_srcdir)/src/libs/libhmat       \
              -I$(top_srcdir)/src/libs/libSGJC       \
              -I$(top_srcdir)/src/libs/libTriInt     \
              -I$(top_srcdir)/src/libs/libhrutil
