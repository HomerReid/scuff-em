EXTRA_DIST = \
 UnitTestSphere.geo	\
 UnitTestSphere.msh  	\
 PECSphere.scuffgeo 	\
 SiSphere.scuffgeo  	\
 TwoSpheres.scuffgeo  

OBJ = $(SRC:.cc=.o)
DAC = $(SRC:.cc=.dac)
DONE = $(SRC:.cc=.done)

.SUFFIXES = .dac .done

LIBSCUFF = $(top_builddir)/src/libs/libscuff/libscuff.la

AM_CPPFLAGS = -DSCUFF \
              -I$(top_srcdir)/src/libs/libscuff      \
              -I$(top_srcdir)/src/libs/libIncField   \
              -I$(top_srcdir)/src/libs/libMatProp    \
              -I$(top_srcdir)/src/libs/libMDInterp    \
              -I$(top_srcdir)/src/libs/libhmat       \
              -I$(top_srcdir)/src/libs/libSGJC       \
              -I$(top_srcdir)/src/libs/libTriInt     \
              -I$(top_srcdir)/src/libs/libhrutil

noinst_PROGRAMS = 		\
 unit-test-compactMatrix	\
 unit-test-ppis

check_PROGRAMS = 		\
 unit-test-compactMatrix	\
 unit-test-ppis

TESTS = 			\
 unit-test-compactMatrix	\
 unit-test-ppis

unit_test_ppis_SOURCES = unit-test-ppis.cc
unit_test_ppis_LDADD = $(LIBSCUFF)

unit_test_compactMatrix_SOURCES = unit-test-compactMatrix.cc
unit_test_compactMatrix_LDADD   = $(LIBSCUFF)
