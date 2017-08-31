noinst_LTLIBRARIES = libSubstrate.la
pkginclude_HEADERS = libSubstrate.h

libSubstrate_la_SOURCES = 	\
 libSubstrate.h		\
 libSubstrate.cc	\
 Static.cc     		

AM_CPPFLAGS = -I$(top_srcdir)/src/libs/libhrutil    \
              -I$(top_srcdir)/src/libs/libhmat      \
              -I$(top_srcdir)/src/libs/libMatProp   \
              -I$(top_srcdir)/src/libs/libMDInterp  \
              -I$(top_srcdir)/src/libs/libSGJC      \
              -I$(top_srcdir)/src/libs/libSpherical \
              -I$(top_builddir) # for config.h

OTHER_LIBS=                      \
 ../libMatProp/libMatProp.la     \
 ../libMDInterp/libMDInterp.la   \
 ../libhmat/libhmat.la           \
 ../libSGJC/libSGJC.la           \
 ../libSpherical/libSpherical.la \
 ../libhrutil/libhrutil.la     

noinst_PROGRAMS = tStaticSubstrate tFullWaveSubstrate

tStaticSubstrate_SOURCES = tStaticSubstrate.cc
tStaticSubstrate_LDADD = libSubstrate.la $(OTHER_LIBS)

tFullWaveSubstrate_SOURCES = tFullWaveSubstrate.cc
tFullWaveSubstrate_LDADD = libSubstrate.la $(OTHER_LIBS)
