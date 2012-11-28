EXTRA_DIST = 					\
 README						\
 UnitTestSphere.geo				\
 UnitTestSphere.msh  				\
 UnitTestSphere_Coarse.msh  			\
 Square.geo        				\
 Square_Coarse.msh				\
 PECSphere.scuffgeo 				\
 SiSphere.scuffgeo 		 		\
 TwoSpheres.scuffgeo 				\
 TwoSpheres_Coarse.scuffgeo			\
 SphereSubstrateArray_Coarse.scuffgeo

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
 unit-test-PBCMatrix     	\
 unit-test-PPIs

check_PROGRAMS = 		\
 unit-test-compactMatrix	\
 unit-test-PBCMatrix     	\
 unit-test-PPIs

TESTS = 			\
 unit-test-compactMatrix	\
 unit-test-PBCMatrix     	\
 unit-test-PPIs

unit_test_PPIs_SOURCES = unit-test-PPIs.cc
unit_test_PPIs_LDADD = $(LIBSCUFF)

unit_test_compactMatrix_SOURCES = unit-test-compactMatrix.cc
unit_test_compactMatrix_LDADD   = $(LIBSCUFF)

unit_test_PBCMatrix_SOURCES = unit-test-PBCMatrix.cc
unit_test_PBCMatrix_LDADD   = $(LIBSCUFF)
