noinst_LTLIBRARIES = libIncField.la
pkginclude_HEADERS = libIncField.h
libIncField_la_SOURCES = IncField.cc GaussianBeam.cc PlaneWave.cc PointSource.cc libIncField.h libVec.h

AM_CPPFLAGS = -I$(top_srcdir)/src/libs/libhrutil

# 20141001 commenting out this test program for now because 
# it depends on libscuff
#noinst_PROGRAMS = tPointSource
#tPointSource_SOURCES = tPointSource.cc
#tPointSource_LDADD = libIncField.la ../libhrutil/libhrutil.la ../libscuff/libscuff.la
