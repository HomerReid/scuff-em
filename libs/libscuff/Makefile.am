lib_LTLIBRARIES = libscuff.la
pkginclude_HEADERS = libscuff.h \
  EquivalentEdgePairs.h		\
  GTransformation.h     	\
  GBarAccelerator.h		\
  PFTOptions.h			\
  PanelCubature.h

libscuff_la_SOURCES =		\
 AssembleBEMMatrix2018.cc      	\
 AssembleBEMMatrix.cc          	\
 AssembleRHSVector.cc 		\
 AssessPanelPair.cc 		\
 CalcGC.cc 			\
 DSIPFT.cc 			\
 EdgeEdgeInteractions.cc	\
 EMTPFT.cc			\
 EPPFT.cc			\
 EquivalentEdgePairs.cc     	\
 EquivalentEdgePairs.h      	\
 ExpandCurrentDistribution.cc 	\
 Faddeeva.cc        		\
 Faddeeva.hh        		\
 FIBBICache.cc   		\
 FIPPICache.cc 			\
 GBarAccelerator.cc 		\
 GBarAccelerator.h  		\
 GBarVDEwald.cc     		\
 GCMatrixElements.cc		\
 GetDipoleMoments.cc 		\
 GetDyadicGFs.cc        	\
 GetFields.cc 			\
 GetNearFields.cc 		\
 GetPFT.cc			\
 GetSphericalMoments.cc 	\
 GTransformation.cc 		\
 GTransformation.h 		\
 InitEdgeList.cc 		\
 libscuff.h 			\
 libscuffInternals.h		\
 MomentPFT.cc			\
 OPFT.cc  			\
 PanelCubature.cc          	\
 PanelCubature.h		\
 PanelPanelInteractions.cc 	\
 PBCSetup.cc 			\
 PFTOptions.h			\
 PointInObject.cc 		\
 QDFIPPI.cc 			\
 QIFIPPICubature.cc 		\
 QIFIPPITaylorDuffy.cc 		\
 QIFIPPITaylorDuffyV2P0.cc 	\
 ReadComsolFile.cc 		\
 ReadGMSHFile.cc 		\
 RWGGeometry.cc 		\
 RWGSurface.cc 			\
 rwlock.cc 			\
 rwlock.h 			\
 SurfaceSurfaceInteractions.cc 	\
 TaylorDuffy.cc 		\
 TaylorDuffy.h 			\
 Visualize.cc 			

# combine all of the auxiliary libraries into a single library 
libscuff_la_LIBADD = 		 \
 ../libhrutil/libhrutil.la	 \
 ../libhmat/libhmat.la		 \
 ../libSGJC/libSGJC.la		 \
 ../libIncField/libIncField.la   \
 ../libMDInterp/libMDInterp.la   \
 ../libMatProp/libMatProp.la	 \
 ../libTriInt/libTriInt.la 	 \
 ../libSpherical/libSpherical.la \
 ../libSubstrate/libSubstrate.la \
 ../libscuffSolver/libscuffSolver.la   \
 ../libStaticSolver/libStaticSolver.la

# set the shared-library version number (DIFFERENT from human version #)
libscuff_la_LDFLAGS = -version-info @SHARED_VERSION_INFO@

AM_CPPFLAGS = -I$(top_srcdir)/libs/libhrutil    \
              -I$(top_srcdir)/libs/libMatProp   \
              -I$(top_srcdir)/libs/libMatProp/cmatheval  \
              -I$(top_srcdir)/libs/libhmat      \
              -I$(top_srcdir)/libs/libSGJC      \
              -I$(top_srcdir)/libs/libIncField  \
              -I$(top_srcdir)/libs/libMDInterp  \
              -I$(top_srcdir)/libs/libSpherical \
              -I$(top_srcdir)/libs/libTriInt    \
              -I$(top_srcdir)/libs/libSubstrate \
              -I$(top_srcdir)/libs/libStaticSolver  \
              -I$(top_srcdir)/libs/libscuffSolver  \
              -I$(top_builddir) # for config.h 
