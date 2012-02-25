AC_INIT(scuff-em, 1.0, homer@homerreid.com)

AC_CONFIG_AUX_DIR([config])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_SRCDIR([README])

AM_INIT_AUTOMAKE

AC_PROG_CC
AC_PROG_CXX
AC_PROG_F77
AC_PROG_RANLIB
AC_PROG_INSTALL

##################################################
# check for pthreads #############################
##################################################
#m4_include([m4/ax_pthread.m4])
#
#AX_PTHREAD([],
#           [AC_MSG_ERROR([pthreads library is required])]
#          )

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
sinclude(acx_blas.m4)
sinclude(acx_lapack.m4)

ACX_BLAS([], [AC_MSG_ERROR([BLAS was not found!])])
ACX_LAPACK([], [AC_MSG_ERROR([LAPACK was not found!])])


LIBS="$LAPACK_LIBS $BLAS_LIBS $LIBS $FLIBS"

#ACX_LAPACK( [],
#            [AC_MSG_ERROR([lapack and blas libraries are required])]
#          )

##################################################
# checks for muparser
##################################################
AC_LANG([C++]) AC_CHECK_HEADER([muParser/muParser.h],[],
                [AC_MSG_ERROR([muParser is required])]
               )
AC_CHECK_LIB( [muparser], 
              [_init], 
              [], 
              [AC_MSG_ERROR([muParser library is required])]
            )

##################################################
##################################################
##################################################
AC_CONFIG_FILES([
 Makefile
 src/Makefile
 src/libs/Makefile
 src/libs/libscuff/Makefile
 src/libs/libIncField/Makefile
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
AC_MSG_NOTICE([** Now type 'make install.'.])
AC_MSG_NOTICE([**])     
