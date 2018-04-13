bin_PROGRAMS = scuff-heat

scuff_heat_SOURCES = 		\
 FrequencyIntegrand.cc 		\
 CreateSHData.cc       		\
 scuff-heat.cc         		\
 scuff-heat.h

scuff_heat_LDADD = $(top_builddir)/libs/libscuff/libscuff.la

AM_CPPFLAGS = -I$(top_srcdir)/libs/libscuff      \
              -I$(top_srcdir)/libs/libIncField   \
              -I$(top_srcdir)/libs/libMatProp    \
              -I$(top_srcdir)/libs/libMDInterp   \
              -I$(top_srcdir)/libs/libhmat       \
              -I$(top_srcdir)/libs/libSGJC       \
              -I$(top_srcdir)/libs/libTriInt     \
              -I$(top_srcdir)/libs/libhrutil
