pkgdatadir = $(datadir)/scuff-em/tests/Fresnel
pkgdata_DATA = E10HalfSpace_40.scuffgeo 		\
               Square_40.msh				\
               E10HalfSpace_40.timing.reference	 	\
               E10HalfSpace_40.transmission.reference	\
               Fresnel.Checklist

AM_TESTS_ENVIRONMENT=\
 export PATH=${PATH}:$(top_builddir)/src/applications/scuff-transmission; \
 export CHECKSCUFFDATA=$(top_buildir)/tests/CheckSCUFFData;               \
 export SCUFF_LOGLEVEL=VERBOSE2;

pkgdata_SCRIPTS = TestFresnel.sh
TESTS = TestFresnel.sh
