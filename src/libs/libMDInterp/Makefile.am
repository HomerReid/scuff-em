noinst_LTLIBRARIES = libMDInterp.la
pkginclude_HEADERS = libMDInterp.h
libMDInterp_la_SOURCES = BinSearch.cc Interp1D.cc Interp2D.cc Interp3D.cc Interp4D.cc freadEC.cc libMDInterp.h

check_PROGRAMS      = tlibMDInterp

tlibMDInterp_SOURCES = tlibMDInterp.cc
tlibMDInterp_LDADD   = ../libscuff/libscuff.la

TESTS = tlibMDInterp

EXTRA_DIST = COPYRIGHT

AM_CPPFLAGS = -I$(top_srcdir)/src/libs/libhrutil	\
              -I$(top_srcdir)/src/libs/libhmat \
              -I$(top_builddir) # for config.h

# bitrotted:  tInterp3D2
#noinst_PROGRAMS = WriteMatrix WriteMatrix2D WriteMatrix56 WriteMatrix64 WriteMatrixLP WriteMatrixLP64 doom tInterp1D tInterp2D tInterp3D
#
#WriteMatrix_SOURCES = WriteMatrix.cc
#WriteMatrix_LDADD = ../libhmat/libhmat.la ../libhrutil/libhrutil.la
#WriteMatrix2D_SOURCES = WriteMatrix2D.cc
#WriteMatrix2D_LDADD = ../libhmat/libhmat.la ../libhrutil/libhrutil.la
#WriteMatrix56_SOURCES = WriteMatrix56.cc
#WriteMatrix56_LDADD = ../libhmat/libhmat.la ../libhrutil/libhrutil.la
#WriteMatrix64_SOURCES = WriteMatrix64.cc
#WriteMatrix64_LDADD = ../libhmat/libhmat.la ../libhrutil/libhrutil.la
#WriteMatrixLP_SOURCES = WriteMatrixLP.cc
#WriteMatrixLP_LDADD = ../libhmat/libhmat.la ../libhrutil/libhrutil.la
#WriteMatrixLP64_SOURCES = WriteMatrixLP64.cc
#WriteMatrixLP64_LDADD = ../libhmat/libhmat.la ../libhrutil/libhrutil.la
#doom_SOURCES = doom.cc
#tInterp1D_SOURCES = tInterp1D.cc
#tInterp1D_LDADD = libMDInterp.la ../libhmat/libhmat.la ../libhrutil/libhrutil.la
#tInterp2D_SOURCES = tInterp2D.cc
#tInterp2D_LDADD = libMDInterp.la ../libhmat/libhmat.la ../libhrutil/libhrutil.la
#tInterp3D_SOURCES = tInterp3D.cc
#tInterp3D_LDADD = libMDInterp.la ../libhmat/libhmat.la ../libhrutil/libhrutil.la
#tInterp3D2_SOURCES = tInterp3D2.cc
#tInterp3D2_LDADD = libMDInterp.la ../libhmat/libhmat.la ../libhrutil/libhrutil.la
