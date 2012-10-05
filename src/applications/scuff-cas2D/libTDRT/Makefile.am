noinst_LTLIBRARIES = libTDRT.la

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
 Visualize.cc			


AM_CPPFLAGS = \
  -I$(top_srcdir)/src/applications/scuff-cas2D/libAmosBessel	\
  -I$(top_srcdir)/src/libs/libSGJC       \
  -I$(top_srcdir)/src/libs/libMDInterp   \
  -I$(top_srcdir)/src/libs/libMatProp    \
  -I$(top_srcdir)/src/libs/libhmat       \
  -I$(top_srcdir)/src/libs/libhrutil	

