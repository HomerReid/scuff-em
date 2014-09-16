lib_LTLIBRARIES = libscuff.la
pkginclude_HEADERS = libscuff.h GTransformation.h FieldGrid.h
libscuff_la_SOURCES = \
 RWGGeometry.cc 		\
 RWGSurface.cc 			\
 ReadComsolFile.cc 		\
 ReadGMSHFile.cc 		\
 InitEdgeList.cc 		\
 PBCSetup.cc 			\
 GTransformation.cc 		\
 GTransformation.h 		\
 PointInObject.cc 		\
 Visualize.cc 			\
 AssembleBEMMatrix.cc          	\
 SurfaceSurfaceInteractions.cc 	\
 EdgeEdgeInteractions.cc	\
 PanelPanelInteractions.cc 	\
 TaylorDuffy.cc 		\
 TaylorDuffy.h 			\
 QDFIPPI.cc 			\
 QIFIPPICubature.cc 		\
 QIFIPPITaylorDuffy.cc 		\
 QIFIPPITaylorDuffyV2P0.cc 	\
 FIPPICache.cc 			\
 GBarAccelerator.cc 		\
 GBarAccelerator.h  		\
 GBarVDEwald.cc     		\
 Faddeeva.cc        		\
 Faddeeva.hh        		\
 FieldGrid.cc 			\
 GetFields.cc 			\
 GetNearFields.cc 		\
 EPPFT.cc 			\
 SIPFT.cc 			\
 OPFT.cc  			\
 GetDipoleMoments.cc 		\
 GetSphericalMoments.cc 	\
 GetDyadicGFs.cc        	\
 ExpandCurrentDistribution.cc 	\
 AssembleRHSVector.cc 		\
 AssessPanelPair.cc 		\
 PanelCubature.cc   		\
 CalcGC.cc 			\
 rwlock.cc 			\
 rwlock.h 			\
 scuffMisc.cc			\
 libscuff.h 			\
 libscuffInternals.h

# bitrotted: AssembleDMDVMatrix.cc BORTObject.cc GetSphericalMoments.cc ParseGeoFiles.cc SphericalMoments.cc

# combine all of the auxiliary libraries into a single library 
libscuff_la_LIBADD = 		\
 ../libhrutil/libhrutil.la	\
 ../libhmat/libhmat.la		\
 ../libSGJC/libSGJC.la		\
 ../libIncField/libIncField.la  \
 ../libMDInterp/libMDInterp.la  \
 ../libMatProp/libMatProp.la	\
 ../libTriInt/libTriInt.la 	\
 ../libSpherical/libSpherical.la

# set the shared-library version number (DIFFERENT from human version #)
libscuff_la_LDFLAGS = -version-info @SHARED_VERSION_INFO@

AM_CPPFLAGS = -I$(top_srcdir)/src/libs/libhrutil    \
              -I$(top_srcdir)/src/libs/libMatProp   \
              -I$(top_srcdir)/src/libs/libMatProp/cmatheval  \
              -I$(top_srcdir)/src/libs/libhmat      \
              -I$(top_srcdir)/src/libs/libSGJC      \
              -I$(top_srcdir)/src/libs/libIncField  \
              -I$(top_srcdir)/src/libs/libMDInterp  \
              -I$(top_srcdir)/src/libs/libSpherical \
              -I$(top_srcdir)/src/libs/libTriInt    \
              -I$(top_builddir) # for config.h 
