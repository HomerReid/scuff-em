lib_LIBRARIES = libscuff.a

include_HEADERS = libscuff.h

libscuff_a_SOURCES =  		\
 AssembleRHSVector.o		\
 AssessPanelPair.o		\
 CalcGC.o			\
 GetDipoleMoments.o		\
 GetFields.o			\
 GTransformation.o		\
 EdgeEdgeInteractions.o     	\
 ExpandCurrentDistribution.o	\
 QDFIPPI.o                    	\
 FIPPICache.o                   \
 QIFIPPICubature.o         	\
 QIFIPPITaylorDuffy.o         	\
 InitEdgeList.o			\
 RWGObject.o			\
 ReadComsolFile.o		\
 ReadGMSHFile.o			\
 PanelPanelInteractions.o	\
 TaylorDuffy.o			\
 RWGGeometry.o			\
 Visualize.o			\
 scuffMisc.o			\
 libscuff.h			\
 libscuffInternals.h		\
 GTransformation.h

AM_CPPFLAGS = -I$(top_srcdir)/src/libs/libhrutil    \
              -I$(top_srcdir)/src/libs/libMatProp   \
              -I$(top_srcdir)/src/libs/libhmat      \
              -I$(top_srcdir)/src/libs/libSGJC      \
              -I$(top_srcdir)/src/libs/libIncField  \
              -I$(top_srcdir)/src/libs/libMDInterp  \
              -I$(top_srcdir)/src/libs/libTriInt
