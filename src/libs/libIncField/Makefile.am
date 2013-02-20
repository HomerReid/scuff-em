noinst_LTLIBRARIES = libIncField.la
pkginclude_HEADERS = libIncField.h
libIncField_la_SOURCES = IncField.cc GaussianBeam.cc PlaneWave.cc PointSource.cc libIncField.h libVec.h

AM_CPPFLAGS = -I$(top_srcdir)/src/libs/libhrutil

noinst_PROGRAMS = tPointSource
tPointSource_SOURCES = tPointSource.cc
tPointSource_LDADD = libIncField.la ../libhrutil/libhrutil.la
