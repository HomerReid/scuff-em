bin_PROGRAMS = scuff-spectrum

scuff_spectrum_SOURCES = scuff-spectrum.cc OutputModules.cc
scuff_spectrum_LDADD = libBeyn.la \
                       $(top_builddir)/src/libs/libscuff/libscuff.la

noinst_LTLIBRARIES = libBeyn.la
libBeyn_la_SOURCES = libBeyn.h libBeyn.cc

noinst_PROGRAMS = tBeyn411
tBeyn411_SOURCES = tBeyn411.cc
tBeyn411_LDADD   = libBeyn.la                                      \
                   $(top_builddir)/src/libs/libhmat/libhmat.la     \
                   $(top_builddir)/src/libs/libhrutil/libhrutil.la

AM_CPPFLAGS = -DSCUFF                                \
              -I$(top_srcdir)/src/libs/libscuff      \
              -I$(top_srcdir)/src/libs/libIncField   \
              -I$(top_srcdir)/src/libs/libMatProp    \
              -I$(top_srcdir)/src/libs/libMDInterp   \
              -I$(top_srcdir)/src/libs/libhmat       \
              -I$(top_srcdir)/src/libs/libSGJC       \
              -I$(top_srcdir)/src/libs/libSubstrate  \
              -I$(top_srcdir)/src/libs/libTriInt     \
              -I$(top_srcdir)/src/libs/libhrutil     \
              -I$(top_builddir) # for config.h
