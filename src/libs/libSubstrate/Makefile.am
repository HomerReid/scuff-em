noinst_LTLIBRARIES = libSubstrate.la
pkginclude_HEADERS = libSubstrate.h

libSubstrate_la_SOURCES = 	\
 libSubstrate.h			\
 libSubstrate.cc		\
 FullWave.cc                	\
 FullWave_SurfaceCurrents.cc	\
 Static.cc     		

# FullWave_PlaneWave.cc		


#AM_CXXFLAGS = -ggdb -O0
#CXXFLAGS = -ggdb -O0
AM_CXXFLAGS = -O3
CXXFLAGS = -O3

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

<<<<<<< HEAD
#noinst_PROGRAMS = tStaticSubstrate tFullWaveSubstrate tGDecomposition
noinst_PROGRAMS = tFullWaveSubstrate
=======
noinst_PROGRAMS = tStaticSubstrate tFullWaveSubstrate 
#tGTwiddle
>>>>>>> 693b027d31d12e8f3090b5360f1ee5b3aac90e75

tStaticSubstrate_SOURCES = tStaticSubstrate.cc
tStaticSubstrate_LDADD = libSubstrate.la $(OTHER_LIBS)

tFullWaveSubstrate_SOURCES = tFullWaveSubstrate.cc
tFullWaveSubstrate_LDADD = libSubstrate.la $(OTHER_LIBS)

<<<<<<< HEAD
tGDecomposition_SOURCES = tGDecomposition.cc
tGDecomposition_LDADD = libSubstrate.la $(OTHER_LIBS)

#tGTwiddle
=======
>>>>>>> 693b027d31d12e8f3090b5360f1ee5b3aac90e75
#tGTwiddle_SOURCES = tGTwiddle.cc
#tGTwiddle_LDADD = libSubstrate.la $(OTHER_LIBS)
