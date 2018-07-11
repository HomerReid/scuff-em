pkgdatadir = $(datadir)/scuff-em/tests/Fresnel
pkgdata_DATA = E10HalfSpace_40.scuffgeo 		\
               Square_40.msh				\
               Fresnel.Checklist

referencedir = $(pkgdatadir)/reference
reference_DATA = reference/E10HalfSpace_40.Runtime      \
                 reference/E10HalfSpace_40.transmission

pkgdata_SCRIPTS = TestFresnel.sh
TESTS = TestFresnel.sh
