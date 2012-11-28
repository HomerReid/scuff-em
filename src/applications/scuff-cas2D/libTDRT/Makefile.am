noinst_LTLIBRARIES = libTDRT.la
noinst_HEADERS = libTDRT.h

libTDRT_la_SOURCES = 			\
 Assemble.cc				\
 HRBesselK.cc				\
 LFunctions.cc				\
 ReadComsolFile.cc			\
 ReadGMSHFile.cc			\
 SCR.cc					\
 StaticSSIDataMap.cc			\
 StaticSSIDataRecord.cc			\
 StaticSSIDataRecordDuffy.cc		\
 TDRTGeometry.cc			\
 TDRTMisc.cc				\
 TDRTObject.cc				\
 uupIntegralCubature.cc			\
 uupIntegralDuffy.cc			\
 uupIntegralSameSegment.cc		\
 Visualize.cc				\
 TValues.cc                  		\
 G1G2Values.cc                  


AM_CPPFLAGS = \
  -I$(top_srcdir)/src/libs/libhrutil	\
  -I$(top_srcdir)/src/applications/scuff-cas2D/libAmosBessel	\
  -I$(top_srcdir)/src/libs/libSGJC       \
  -I$(top_srcdir)/src/libs/libMDInterp   \
  -I$(top_srcdir)/src/libs/libMatProp    \
  -I$(top_srcdir)/src/libs/libhmat

