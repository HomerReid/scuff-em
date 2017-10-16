noinst_LTLIBRARIES = libSubstrate.la
pkginclude_HEADERS = libSubstrate.h

libSubstrate_la_SOURCES = 	\
 libSubstrate.h			\
 libSubstrate.cc		\
 FullWave.cc                	\
 FullWave_SurfaceCurrents.cc	\
 Static.cc     		

#if $(MMAKE)
#AM_CXXFLAGS=-ggdb -O0
#AM_CFLAGS=-ggdb -O0
#  echo "Building with debug support."
#else
  AM_CXXFLAGS = -O3
#  echo "Building without debug support."
#endif

AM_CPPFLAGS = -I$(top_srcdir)/src/libs/libhrutil    \
              -I$(top_srcdir)/src/libs/libhmat      \
              -I$(top_srcdir)/src/libs/libMatProp   \
              -I$(top_srcdir)/src/libs/libMDInterp  \
              -I$(top_srcdir)/src/libs/libSGJC      \
              -I$(top_srcdir)/src/libs/libSpherical \
              -I$(top_srcdir)/src/libs/libIncField  \
              -I$(top_srcdir)/src/libs/libscuff     \
              -I$(top_builddir) # for config.h

OTHER_LIBS=                      \
 ../libMatProp/libMatProp.la     \
 ../libMDInterp/libMDInterp.la   \
 ../libhmat/libhmat.la           \
 ../libSGJC/libSGJC.la           \
 ../libSpherical/libSpherical.la \
 ../libhrutil/libhrutil.la       \
 ../libscuff/libscuff.la      

noinst_PROGRAMS = tStaticSubstrate tFullWaveSubstrate
check_PROGRAMS = tStaticSubstrate tFullWaveSubstrate tGTwiddle tRotateG

tStaticSubstrate_SOURCES = tStaticSubstrate.cc
tStaticSubstrate_LDADD = libSubstrate.la $(OTHER_LIBS)

tFullWaveSubstrate_SOURCES = tFullWaveSubstrate.cc
tFullWaveSubstrate_LDADD = libSubstrate.la $(OTHER_LIBS)

tGTwiddle_SOURCES = tGTwiddle.cc
tGTwiddle_LDADD = libSubstrate.la $(OTHER_LIBS)

tRotateG_SOURCES = tRotateG.cc
tRotateG_LDADD = libSubstrate.la $(OTHER_LIBS)
