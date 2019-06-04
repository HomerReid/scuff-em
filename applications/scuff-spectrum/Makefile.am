bin_PROGRAMS = scuff-spectrum

scuff_spectrum_SOURCES = scuff-spectrum.cc OutputModules.cc
scuff_spectrum_LDADD = libBeyn.la \
                       $(top_builddir)/libs/libscuff/libscuff.la

noinst_LTLIBRARIES = libBeyn.la
libBeyn_la_SOURCES = libBeyn.h libBeyn.cc

noinst_PROGRAMS = tBeyn411
tBeyn411_SOURCES = tBeyn411.cc
tBeyn411_LDADD   = libBeyn.la                                      \
                   $(top_builddir)/libs/libhmat/libhmat.la     \
                   $(top_builddir)/libs/libhrutil/libhrutil.la

AM_CPPFLAGS = -DSCUFF                                \
              -DVERSIONSTRING=\"$(shell git describe always --dirty)\"\
              -I$(top_srcdir)/libs/libscuff      \
              -I$(top_srcdir)/libs/libIncField   \
              -I$(top_srcdir)/libs/libMatProp    \
              -I$(top_srcdir)/libs/libMDInterp   \
              -I$(top_srcdir)/libs/libhmat       \
              -I$(top_srcdir)/libs/libSGJC       \
              -I$(top_srcdir)/libs/libSubstrate  \
              -I$(top_srcdir)/libs/libTriInt     \
              -I$(top_srcdir)/libs/libhrutil     \
              -I$(top_builddir) # for config.h    
