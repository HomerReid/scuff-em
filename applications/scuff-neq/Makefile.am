bin_PROGRAMS = scuff-neq scuff-integrate

scuff_neq_SOURCES = 		\
 scuff-neq.cc         		\
 scuff-neq.h			\
 WriteFlux.cc             	\
 CreateSNEQData.cc

scuff_neq_LDADD = $(top_builddir)/libs/libscuff/libscuff.la

scuff_integrate_SOURCES = scuff-integrate.cc

scuff_integrate_LDADD = $(top_builddir)/libs/libscuff/libscuff.la

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
