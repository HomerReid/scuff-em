lib_LTLIBRARIES = libscuff.la
pkginclude_HEADERS = libscuff.h GTransformation.h
libscuff_la_SOURCES = AssembleBEMMatrix.cc AssembleRHSVector.cc AssessPanelPair.cc CalcGC.cc EdgeEdgeInteractions.cc ExpandCurrentDistribution.cc FIPPICache.cc GTransformation.cc GTransformation.h GetDipoleMoments.cc GetFields.cc InitEdgeList.cc PanelPanelInteractions.cc QDFIPPI.cc QIFIPPICubature.cc QIFIPPITaylorDuffy.cc RWGGeometry.cc RWGObject.cc ReadComsolFile.cc ReadGMSHFile.cc TaylorDuffy.cc TaylorDuffy.h Visualize.cc libscuff.h libscuffInternals.h scuffMisc.cc

# bitrotted: AssembleDMDVMatrix.cc BORTObject.cc GetSphericalMoments.cc ParseGeoFiles.cc SphericalMoments.cc

# combine all of the auxiliary libraries into a single library
libscuff_la_LIBADD = ../libhrutil/libhrutil.la ../libhmat/libhmat.la	\
../libSGJC/libSGJC.la ../libIncField/libIncField.la			\
../libMDInterp/libMDInterp.la ../libMatProp/libMatProp.la		\
../libTriInt/libTriInt.la

AM_CPPFLAGS = -I$(top_srcdir)/src/libs/libhrutil    \
              -I$(top_srcdir)/src/libs/libMatProp   \
              -I$(top_srcdir)/src/libs/libhmat      \
              -I$(top_srcdir)/src/libs/libSGJC      \
              -I$(top_srcdir)/src/libs/libIncField  \
              -I$(top_srcdir)/src/libs/libMDInterp  \
              -I$(top_srcdir)/src/libs/libTriInt
