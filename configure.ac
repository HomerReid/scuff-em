AC_INIT(scuff-em, 1.0, homer@homerreid.com)
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

AC_PROG_CC
AC_PROG_CXX
AC_PROG_F77
AC_PROG_INSTALL

##################################################
# checks for hdf5
##################################################
AC_CHECK_LIB( [hdf5],
              [H5Fcreate],
              [],
              [AC_MSG_ERROR([hdf5 library is required])]
            )

AC_CHECK_LIB( [hdf5_hl],
              [H5LTmake_dataset_double],
              [],
              [AC_MSG_ERROR([hdf5 high-level library is required])]
            )

##################################################
# checks for blas/lapack
##################################################

ACX_BLAS([], [AC_MSG_ERROR([BLAS was not found!])])
ACX_LAPACK([], [AC_MSG_ERROR([LAPACK was not found!])])

LIBS="$LAPACK_LIBS $BLAS_LIBS $LIBS $FLIBS"

##################################################
# checks for muparser
##################################################
AC_LANG([C++])
AC_CHECK_HEADERS([muParser/muParser.h, muParser.h])
save_LIBS=$LIBS
LIBS="$LIBS -lmuparser"
AC_MSG_CHECKING([for -lmuparser])
AC_LINK_IFELSE([AC_LANG_PROGRAM([
#ifdef HAVE_MUPARSER_H
#  include <muParser.h>
#else
#  include <muParser/muParser.h>
#endif], [mu::Parser p;])],
               [AC_DEFINE([HAVE_MUPARSER],[1],[have muParser library])
                AC_MSG_RESULT([yes])],
               [AC_MSG_RESULT([no])
                AC_MSG_WARN([muParser library not found, user-defined material functions will not be supported])
                 LIBS=$save_LIBS]
              )

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
 src/libs/libSGJC/Makefile
 src/libs/libTriInt/Makefile
 src/libs/libhmat/Makefile
 src/libs/libhrutil/Makefile
 src/libs/libscuff/Makefile
 src/applications/Makefile
 src/applications/scuff-scatter/Makefile
 src/applications/scuff-caspol/Makefile
 src/applications/scuff-heat/Makefile
 src/applications/scuff-analyze/Makefile
])
AC_OUTPUT

##################################################
##################################################
##################################################
AC_MSG_NOTICE([**])
AC_MSG_NOTICE([** The configure process was successful.])
AC_MSG_NOTICE([** Now type 'make' to compile and sudo 'make install'.])
AC_MSG_NOTICE([**])     
