AC_INIT(scuff-em, 0.95, homer@homerreid.com)
AC_CONFIG_SRCDIR([README])

# Shared-library version.  This has a technical
# meaning related to binary compatibility and is NOT
# the same as the human version number.  See http://tinyurl.com/yegqrlz
SHARED_VERSION_INFO="0:0:0" # CURRENT:REVISION:AGE

AM_INIT_AUTOMAKE
LT_INIT
AM_CONFIG_HEADER([config.h])
AC_CONFIG_MACRO_DIR([m4])
AM_MAINTAINER_MODE
AC_SUBST(SHARED_VERSION_INFO)
AC_DISABLE_SHARED dnl shared libs by default cause too much trouble

AC_PROG_CXX
AC_PROG_F77
AC_F77_WRAPPERS
AC_PROG_INSTALL

AC_PROG_YACC
AM_PROG_LEX

AC_LANG([C++])

##################################################
# compiler flags
##################################################

AC_ARG_ENABLE(debug, [AC_HELP_STRING([--enable-debug],[compile for debugging])], enable_debug=$enableval, enable_debug=no)

if test "x$enable_debug" = xyes; then
   CXXFLAGS=-g
   FFLAGS=-g
fi

# TODO: set optimization flags? Does it matter, or is the default -g -O2 okay?

# add gcc warnings in debug or maintainer mode
if test "x$enable_debug" = xyes || test "$USE_MAINTAINER_MODE" = yes; then
if test "$ac_test_CFLAGS" != "set"; then
   if test $ac_cv_prog_gcc = yes; then
      CXXFLAGS="$CXXFLAGS -Wall -W -Wcast-qual -Wpointer-arith -Wcast-align -pedantic -Wno-long-long -Wshadow -Wwrite-strings -Wredundant-decls -pedantic"
   fi
fi
fi

##################################################
# threading library
##################################################

AC_ARG_WITH(openmp, [AC_HELP_STRING([--without-openmp],[do not use OpenMP directives for parallelism])], with_openmp=$withval, with_openmp=yesdefault)

AC_ARG_WITH(pthreads, [AC_HELP_STRING([--with-pthreads],[use POSIX threads for parallelism])], with_pthreads=$withval, with_pthreads=no)

need_pthreads=no

if test "x$with_pthreads" = xyes; then
   AC_DEFINE([USE_PTHREAD], [1], [Define to use pthread threading.])
   if test "x$with_openmp" = xyes; then
      AC_MSG_ERROR([--with-openmp and --with-pthreads are incompatible])
   fi
   with_openmp=no
   need_pthreads=yes
fi

if test "x$with_openmp" = xyesdefault; then with_openmp=yes; fi
if test "x$with_openmp" = xyes; then
   AC_DEFINE([USE_OPENMP], [1], [Define to use OpenMP threading.])
   AX_OPENMP([], [AC_MSG_ERROR([Could not find OpenMP flag; configure with --without-openmp and/or --with-pthreads])])
   CXXFLAGS="$CXXFLAGS $OPENMP_CXXFLAGS"
   need_pthreads=yes # for pthread_rwlock
fi

AC_MSG_CHECKING([suppatage foryaf 3])

if test $need_pthreads = yes; then
   AX_PTHREAD([], [AC_MSG_ERROR([Could not find pthreads library.])])
   LIBS="$PTHREAD_LIBS $LIBS"
   CXXFLAGS="$CXXFLAGS $PTHREAD_CFLAGS"
fi

##################################################
# checks for hdf5
##################################################

AC_ARG_WITH(hdf5, [AC_HELP_STRING([--with-hdf5=<dir>/--without-hdf5],[look for HDF5 in <dir>/include and <dir>/lib])], with_hdf5=$withval, with_hdf5=yes)

if test "x$with_hdf5" != xno; then
   if test -d "$with_hdf5"; then
     CPPFLAGS="$CPPFLAGS -I$with_hdf5/include"
     LDFLAGS="$LDFLAGS -L$with_hdf5/lib"
   fi

   have_hdf5=no
   AC_CHECK_HEADER([hdf5.h],
     [AC_CHECK_LIB([hdf5], [H5Fcreate],
       [LIBS="-lhdf5 $LIBS";
        AC_CHECK_LIB([hdf5_hl], [H5LTmake_dataset_double], 
	             [LIBS="-lhdf5_hl $LIBS"; have_hdf5=yes])])],
     [AC_CHECK_HEADER([H5pubconf.h],
       [AC_PREPROC_IFELSE([AC_LANG_PROGRAM([#include <H5pubconf.h>
#ifdef H5_HAVE_PARALLEL
#  include <mpi.h>
#endif], [])], [], [AC_MSG_WARN([HDF5 compiled with parallel support requires MPI compiler, e.g. configure with CXX=mpic++])])])])
   if test $have_hdf5 = yes; then
      AC_DEFINE([HAVE_HDF5],[1],[define if we have HDF5])
      LIBS="$LIBS -lhdf5_hl -lhdf5"
   else
      AC_MSG_ERROR([couldn't find HDF5; configure --without-hdf5])
   fi
fi

##################################################
# checks for blas/lapack
##################################################

ACX_BLAS([], [AC_MSG_ERROR([BLAS was not found!])])
ACX_LAPACK([], [AC_MSG_ERROR([LAPACK was not found!])])

LIBS="$LAPACK_LIBS $BLAS_LIBS $LIBS $FLIBS"

##################################################
# checks for readline 
# (which is used by some test programs)
##################################################
AC_CHECK_LIB(readline, readline)

##################################################
# checks for backtrace functionality, used to 
# output useful debugging info to log files in
# the event of code crashes
##################################################
AC_CHECK_HEADERS([execinfo.h])
AC_CHECK_FUNCS([backtrace])

##################################################
# Python SWIG wrappers
##################################################

AC_ARG_WITH(python,
	[AC_HELP_STRING([--without-python], [don't compile Python interface])],
        with_python=$withval,with_python=yes)

if test "$enable_shared" = no; then
   AC_MSG_WARN([Python interface requires --enable-shared; disabling])
   have_python=no
else

if test "x$with_python" = xno; then
  have_python=no
else

AM_PATH_PYTHON([],[have_python=yes],[have_python=no])
if test $have_python = yes; then
  AC_MSG_CHECKING([for Python include directory])
  pinc=`echo "import distutils.sysconfig; print distutils.sysconfig.get_python_inc()" | $PYTHON - 2>/dev/null`
  AC_MSG_RESULT([${pinc:-unknown}])
  test -n "$pinc" && PYTHON_INCLUDES="-I$pinc"
  save_CPPFLAGS=$CPPFLAGS
  CPPFLAGS="$CPPFLAGS $PYTHON_INCLUDES"
  AC_CHECK_HEADER([Python.h], [], [AC_MSG_WARN([disabling Python wrappers]) 
                                   have_python=no])
  if test $have_python = yes; then
    AC_MSG_CHECKING([for Numpy include directory])
    pinc=`echo "import numpy; print numpy.get_include()" | $PYTHON - 2>/dev/null`
    AC_MSG_RESULT([${pinc:-unknown}])
    test -n "$pinc" && PYTHON_INCLUDES="$PYTHON_INCLUDES -I$pinc"
    CPPFLAGS="$save_CPPFLAGS $PYTHON_INCLUDES"
    AC_CHECK_HEADER([numpy/arrayobject.h],[],[
      AC_MSG_WARN([disabling Python wrappers])
      have_python=no],[#include <Python.h>])
  fi
  CPPFLAGS=$save_CPPFLAGS
fi

fi # with_python

fi # if enable_shared

AC_SUBST(PYTHON_INCLUDES)
AM_CONDITIONAL(WITH_PYTHON, test x"$have_python" = "xyes")

##################################################
##################################################
##################################################

AC_CONFIG_FILES([
 Makefile
 m4/Makefile
 src/Makefile
 src/libs/Makefile
 src/libs/libIncField/Makefile
 src/libs/libMDInterp/Makefile
 src/libs/libMatProp/Makefile
 src/libs/libMatProp/cmatheval/Makefile
 src/libs/libSGJC/Makefile
 src/libs/libTriInt/Makefile
 src/libs/libhmat/Makefile
 src/libs/libhrutil/Makefile
 src/libs/libscuff/Makefile
 src/libs/python/Makefile
 src/applications/Makefile
 src/applications/scuff-analyze/Makefile
 src/applications/scuff-caspol/Makefile
 src/applications/scuff-cas2D/libTDRT/Makefile
 src/applications/scuff-cas2D/libAmosBessel/Makefile
 src/applications/scuff-cas2D/Makefile
 src/applications/scuff-cas3D/Makefile
 src/applications/scuff-heat/Makefile
 src/applications/scuff-neq/Makefile
 src/applications/scuff-plotEpsMu/Makefile
 src/applications/scuff-rf/Makefile
 src/applications/scuff-scatter/Makefile
 src/applications/scuff-transmission/Makefile
 examples/Makefile
 examples/ApertureDiffraction/Makefile
 examples/GetTMatrix/Makefile
 examples/PlateWithHole/Makefile
 examples/PerforatedThinFilm/Makefile
 examples/SolidSphere/Makefile
 examples/SphericalShell/Makefile
 examples/TwoSpheres/Makefile
 examples/TwoDimensionalCylinders/Makefile
 examples/ThinFilm/Makefile
 examples/WireAntenna/Makefile
 interactiveTests/Makefile
 unitTests/Makefile
])
AC_OUTPUT

##################################################
##################################################
##################################################
AC_MSG_NOTICE([**])
AC_MSG_NOTICE([** The configure process was successful.])
AC_MSG_NOTICE([** Now type 'make' to compile and sudo 'make install'.])
AC_MSG_NOTICE([**])     
