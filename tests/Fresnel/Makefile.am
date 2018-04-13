pkgdatadir = $(datadir)/scuff-em/tests/Fresnel
pkgdata_DATA = E10HalfSpace_40.scuffgeo 		\
               Square_40.msh				\
               E10HalfSpace_40.timing.reference	 	\
               E10HalfSpace_40.transmission.reference	\
               Fresnel.Checklist

AM_TESTS_ENVIRONMENT=export SCUFFBIN=$(bindir); export SCUFFDATA=$(datadir);

pkgdata_SCRIPTS = TestFresnel.sh
TESTS = TestFresnel.sh
