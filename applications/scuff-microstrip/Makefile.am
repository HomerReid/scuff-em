bin_PROGRAMS = scuff-microstrip

scuff_microstrip_SOURCES =	\
 scuff-microstrip.cc		

scuff_microstrip_LDADD = $(top_builddir)/libs/libscuff/libscuff.la

AM_CPPFLAGS = -DSCUFF \
              -I$(top_srcdir)/libs/libhrutil	\
              -I$(top_srcdir)/libs/libscuff		 \
              -I$(top_srcdir)/libs/libIncField   \
              -I$(top_srcdir)/libs/libMatProp    \
              -I$(top_srcdir)/libs/libMDInterp    \
              -I$(top_srcdir)/libs/libhmat       \
              -I$(top_srcdir)/libs/libSGJC       \
              -I$(top_srcdir)/libs/libSubstrate  \
              -I$(top_srcdir)/libs/libscuffSolver \
              -I$(top_srcdir)/libs/libTriInt
