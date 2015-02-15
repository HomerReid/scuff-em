EXTRA_DIST = 					\
 Sphere_501.msh					\
 PECSphere_501.scuffgeo				\
 SiO2Sphere_501.scuffgeo			\
 GoldSphere_501.scuffgeo			\
 SSphere_255.msh				\
 SSphere_R0P25_255.msh				\
 UnitTestSphere_R0P75_414.msh     		\
 Square_40.msh                    		\
 PECSphere_255.scuffgeo				\
 PECSpheres_255.scuffgeo			\
 SiSphere_255.scuffgeo				\
 SiSpheres_255.scuffgeo				\
 PECSphere_R0P75_414.scuffgeo			\
 PECPlate_40.scuffgeo             		\
 SiSlab_40.scuffgeo               		\
 SphereSlabArray.scuffgeo

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
 unit-test-BEMMatrix     	\
 unit-test-PPIs			\
 unit-test-PFT 

check_PROGRAMS = 		\
 unit-test-BEMMatrix     	\
 unit-test-PPIs			\
 unit-test-PFT

TESTS = 			\
 unit-test-BEMMatrix     	\
 unit-test-PPIs			\
 unit-test-PFT

unit_test_BEMMatrix_SOURCES = unit-test-BEMMatrix.cc
unit_test_BEMMatrix_LDADD   = $(LIBSCUFF)

unit_test_PPIs_SOURCES = unit-test-PPIs.cc
unit_test_PPIs_LDADD = $(LIBSCUFF)

unit_test_PFT_SOURCES = unit-test-PFT.cc
unit_test_PFT_LDADD = $(LIBSCUFF)
