AM_CPPFLAGS = -I$(top_srcdir)/src/libs/libhrutil    \
              -I$(top_srcdir)/src/libs/libhmat      \
              -I$(top_srcdir)/src/libs/libMatProp   \
              -I$(top_srcdir)/src/libs/libMDInterp  \
              -I$(top_srcdir)/src/libs/libSGJC      \
              -I$(top_srcdir)/src/libs/libSpherical \
              -I$(top_srcdir)/src/libs/libIncField  \
              -I$(top_srcdir)/src/libs/libSubstrate \
              -I$(top_srcdir)/src/libs/libscuff     \
              -I$(top_builddir) # for config.h

check_PROGRAMS = tGTwiddle \
 tFullWaveSubstrate        \
 tVTwiddle                 \
 tMOISGFs 
 
TESTS = tGTwiddle tFullWaveSubstrate

TESTLIBS = ../libSubstrate.la      \
../../libscuff/libscuff.la         \
../../libMatProp/libMatProp.la     \
../../libMDInterp/libMDInterp.la   \
../../libhmat/libhmat.la           \
../../libSGJC/libSGJC.la           \
../../libSpherical/libSpherical.la \
../../libhrutil/libhrutil.la

tMOISGFs_SOURCES = tMOISGFs.cc
tMOISGFs_LDADD = $(TESTLIBS)

tVTwiddle_SOURCES = tVTwiddle.cc
tVTwiddle_LDADD = $(TESTLIBS)

tGTwiddle_SOURCES = tGTwiddle.cc
tGTwiddle_LDADD = $(TESTLIBS)
 
tFullWaveSubstrate_SOURCES = tFullWaveSubstrate.cc
tFullWaveSubstrate_LDADD = $(TESTLIBS)

tDerivatives_SOURCES = tDerivatives.cc
tDerivatives_LDADD = $(TESTLIBS)

tqIntegrandPPI_SOURCES = tqIntegrandPPI.cc
tqIntegrandPPI_LDADD = $(TESTLIBS)

tqIntegrandFields_SOURCES = tqIntegrandFields.cc
tqIntegrandFields_LDADD = $(TESTLIBS)
