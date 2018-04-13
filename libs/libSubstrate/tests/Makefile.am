AM_CPPFLAGS = -I$(top_srcdir)/libs/libhrutil    \
              -I$(top_srcdir)/libs/libhmat      \
              -I$(top_srcdir)/libs/libMatProp   \
              -I$(top_srcdir)/libs/libMDInterp  \
              -I$(top_srcdir)/libs/libSGJC      \
              -I$(top_srcdir)/libs/libSpherical \
              -I$(top_srcdir)/libs/libIncField  \
              -I$(top_srcdir)/libs/libSubstrate \
              -I$(top_srcdir)/libs/libscuff     \
              -I$(top_builddir) # for config.h

#check_PROGRAMS = tGTwiddle \
# tFullWaveSubstrate        \
# tMOISGFs                  \
# tInterpolation
 
#TESTS = tGTwiddle tFullWaveSubstrate
TESTS =

TESTLIBS = ../libSubstrate.la      \
../../libscuff/libscuff.la         \
../../libMatProp/libMatProp.la     \
../../libMDInterp/libMDInterp.la   \
../../libhmat/libhmat.la           \
../../libSGJC/libSGJC.la           \
../../libSpherical/libSpherical.la \
../../libhrutil/libhrutil.la
#
#tMOISGFs_SOURCES = tMOISGFs.cc
#tMOISGFs_LDADD = $(TESTLIBS)
#
#tInterpolation_SOURCES = tInterpolation.cc
#tInterpolation_LDADD = $(TESTLIBS)
#
#tVTwiddle_SOURCES = tVTwiddle.cc
#tVTwiddle_LDADD = $(TESTLIBS)
#
#tGTwiddle_SOURCES = tGTwiddle.cc
#tGTwiddle_LDADD = $(TESTLIBS)
 
#tFullWaveSubstrate_SOURCES = tFullWaveSubstrate.cc
#tFullWaveSubstrate_LDADD = $(TESTLIBS)

#tDerivatives_SOURCES = tDerivatives.cc
#tDerivatives_LDADD = $(TESTLIBS)
#
#tqIntegrandPPI_SOURCES = tqIntegrandPPI.cc
#tqIntegrandPPI_LDADD = $(TESTLIBS)
#
#tqIntegrandFields_SOURCES = tqIntegrandFields.cc
#tqIntegrandFields_LDADD = $(TESTLIBS)
