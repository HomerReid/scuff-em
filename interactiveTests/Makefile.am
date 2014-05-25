bin_PROGRAMS = 		\
 scuff-test-overlap
# scuff-test-Ewald   

scuff_test_overlap_SOURCES =	\
 readlineReplacement.h          \
 scuff-test-overlap.cc

#scuff_test_Ewald_SOURCES =	\
# readlineReplacement.h          \
# scuff-test-Ewald.cc

scuff_test_overlap_LDADD     = $(top_builddir)/src/libs/libscuff/libscuff.la

#scuff_test_Ewald_LDADD     = $(top_builddir)/src/libs/libscuff/libscuff.la

AM_CPPFLAGS = -DSCUFF \
              -I$(top_srcdir)/src/libs/libscuff      \
              -I$(top_srcdir)/src/libs/libIncField   \
              -I$(top_srcdir)/src/libs/libMatProp    \
              -I$(top_srcdir)/src/libs/libMDInterp    \
              -I$(top_srcdir)/src/libs/libhmat       \
              -I$(top_srcdir)/src/libs/libSGJC       \
              -I$(top_srcdir)/src/libs/libTriInt     \
              -I$(top_srcdir)/src/libs/libhrutil
