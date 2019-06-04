bin_PROGRAMS = scuff-scatter

scuff_scatter_SOURCES = 	\
 scuff-scatter.cc 		\
 OutputModules.cc 		\
 scuff-scatter.h

scuff_scatter_LDADD = $(top_builddir)/libs/libscuff/libscuff.la

AM_CPPFLAGS = -DSCUFF \
              -I$(top_srcdir)/libs/libscuff      \
              -I$(top_srcdir)/libs/libIncField   \
              -I$(top_srcdir)/libs/libMatProp    \
              -I$(top_srcdir)/libs/libMDInterp    \
              -I$(top_srcdir)/libs/libhmat       \
              -I$(top_srcdir)/libs/libSGJC       \
              -I$(top_srcdir)/libs/libSubstrate  \
              -I$(top_srcdir)/libs/libTriInt     \
              -I$(top_srcdir)/libs/libhrutil     \
              -DVERSIONSTRING=\"$(shell git describe always --dirty)\"
