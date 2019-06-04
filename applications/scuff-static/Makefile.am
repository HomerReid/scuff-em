bin_PROGRAMS = scuff-static

scuff_static_SOURCES = 		\
 scuff-static.h			\
 BMAccelerator.cc		\
 OutputModules.cc		\
 scuff-static.cc

scuff_static_LDADD = $(top_builddir)/libs/libscuff/libscuff.la

AM_CPPFLAGS = -I$(top_srcdir)/libs/libscuff      \
              -I$(top_srcdir)/libs/libSubstrate  \
              -I$(top_srcdir)/libs/libIncField   \
              -I$(top_srcdir)/libs/libMatProp    \
              -I$(top_srcdir)/libs/libMatProp/cmatheval \
              -I$(top_srcdir)/libs/libMDInterp   \
              -I$(top_srcdir)/libs/libhmat       \
              -I$(top_srcdir)/libs/libSGJC       \
              -I$(top_srcdir)/libs/libSpherical  \
              -I$(top_srcdir)/libs/libSubstrate  \
              -I$(top_srcdir)/libs/libTriInt     \
              -I$(top_srcdir)/libs/libStaticSolver \
              -I$(top_srcdir)/libs/libhrutil     \
              -DVERSIONSTRING=\"$(shell git describe always --dirty)\"
