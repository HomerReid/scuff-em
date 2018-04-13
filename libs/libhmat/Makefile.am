noinst_LTLIBRARIES = libhmat.la
pkginclude_HEADERS = libhmat.h

libhmat_la_SOURCES = 	\
 lapack.h		\
 lapack_names.h		\
 LBWrappers.cc 		\
 C2ML.cc 		\
 HDF5IO.cc 		\
 GetEntries.cc		\
 HMatrix.cc 		\
 HVector.cc 		\
 SMatrix.cc		\
 Sort.cc 		\
 TextIO.cc

AM_CPPFLAGS = -I$(top_srcdir)/libs/libhrutil \
              -I$(top_builddir) # for config.h

AM_CXXFLAGS = -Wsign-compare

# 20160522 removing tInvert from the noinst_PROGRAMS because
# it is misbehavin' due to my recent switchover from drand48() to 
# the modern glibc random-number generator, which apparently 
# is not standard across all platforms :(
# tInvert_SOURCES = tInvert.cc
# tInvert_LDADD = libhmat.la ../libhrutil/libhrutil.la

noinst_PROGRAMS = tLUSolve tMultiply tReadFromFile tTextIO tlibhmat2 tQR tGetEntries tSMatrix
tQR_SOURCES = tQR.cc
tQR_LDADD = libhmat.la ../libhrutil/libhrutil.la
tLUSolve_SOURCES = tLUSolve.cc
tLUSolve_LDADD = libhmat.la ../libhrutil/libhrutil.la
tMultiply_SOURCES = tMultiply.cc
tMultiply_LDADD = libhmat.la ../libhrutil/libhrutil.la
tReadFromFile_SOURCES = tReadFromFile.cc
tReadFromFile_LDADD = libhmat.la ../libhrutil/libhrutil.la
tTextIO_SOURCES = tTextIO.cc
tTextIO_LDADD = libhmat.la ../libhrutil/libhrutil.la
tlibhmat2_SOURCES = tlibhmat2.cc
tlibhmat2_LDADD = libhmat.la ../libhrutil/libhrutil.la
tGetEntries_SOURCES = tGetEntries.cc
tGetEntries_LDADD = libhmat.la ../libhrutil/libhrutil.la
tSMatrix_SOURCES = tSMatrix.cc
tSMatrix_LDADD = libhmat.la ../libhrutil/libhrutil.la

BUILT_SOURCES = lapack_names.h

if MAINTAINER_MODE

lapack_names.h: lapack.h
	(echo "/* AUTOMATICALLY GENERATED - DO NOT EDIT */"; echo "#ifndef LAPACK_NAMES_H"; echo "#define LAPACK_NAMES_H 1"; grep '_(' lapack.h |perl -pe 's/.* ([a-z0-9]*)_\(.*$$/\1/' |while read n; do echo "#define $${n}_ F77_FUNC($$n,"`echo $$n | tr a-z A-Z`")"; done; echo "#endif") > $@

endif
