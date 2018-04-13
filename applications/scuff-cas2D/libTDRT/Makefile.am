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
  -I$(top_srcdir)/libs/libhrutil	\
  -I$(top_srcdir)/applications/scuff-cas2D/libAmosBessel	\
  -I$(top_srcdir)/libs/libSGJC       \
  -I$(top_srcdir)/libs/libMDInterp   \
  -I$(top_srcdir)/libs/libMatProp    \
  -I$(top_srcdir)/libs/libhmat

