##################################################
# We build the binary code 'CheckSCUFFData' and
# install it at the top of the test-suite folder
# hierarchy (i.e. $PREFIX/share/scuff-em/tests)
# so it will be available for users to run
# individual tests
##################################################
my_pkgdatadir = $(datadir)/scuff-em/tests
my_pkgdata_PROGRAMS = CheckSCUFFData

CheckSCUFFData_SOURCES = CheckSCUFFData.cc
CheckSCUFFData_LDADD   = $(top_builddir)/libs/libscuff/libscuff.la

AM_CPPFLAGS = -DSCUFF 				\
              -I$(top_srcdir)/libs/libhmat      \
              -I$(top_srcdir)/libs/libhrutil

##################################################
# The makefile for each subdirectory just copies
# files needed for the test to the appropriate
# subdirectory of $PREFIX/share/scuff-em/tests
##################################################
SUBDIRS = Mie Fresnel
