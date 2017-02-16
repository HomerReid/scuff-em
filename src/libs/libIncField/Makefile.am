noinst_LTLIBRARIES = libIncField.la
pkginclude_HEADERS = libIncField.h libVec.h
libIncField_la_SOURCES = \
 IncField.cc      \
 IncFieldList.cc  \
 GaussianBeam.cc  \
 PlaneWave.cc     \
 PointSource.cc   \
 SphericalWave.cc \
 libIncField.h    \
 libVec.h 

AM_CPPFLAGS = -I$(top_srcdir)/src/libs/libhrutil    \
              -I$(top_srcdir)/src/libs/libTriInt    \
              -I$(top_srcdir)/src/libs/libSpherical \
              -I$(top_srcdir)/src/libs/libSGJC      \
              -I$(top_srcdir)/src/libs/libhmat

# 20141001 commenting out this test program for now because 
# it depends on libscuff
#noinst_PROGRAMS = tPointSource
#tPointSource_SOURCES = tPointSource.cc
#tPointSource_LDADD = libIncField.la ../libhrutil/libhrutil.la ../libscuff/libscuff.la
