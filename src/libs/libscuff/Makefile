######################################################################
# include system-dependent paths etc.
#######################################################################
include $(HOME)/work/Makefile.Rules

##################################################
##################################################
##################################################
HEADERS = libscuff.h			\
          libscuffInternals.h		\
          GTransformation.h

OBJS = AssembleBEMMatrix.o		\
       AssembleRHSVector.o		\
       AssessPanelPair.o		\
       CalcGC.o				\
       GetDipoleMoments.o		\
       GetSphericalMoments.o		\
       GetFields.o			\
       GTransformation.o		\
       EdgeEdgeInteractions.o     	\
       ExpandCurrentDistribution.o	\
       QDFIPPI.o                    	\
       FIPPITable.o                    	\
       QIFIPPICubature.o         	\
       QIFIPPITaylorDuffy.o         	\
       InitEdgeList.o			\
       RWGObject.o			\
       ReadComsolFile.o			\
       ReadGMSHFile.o			\
       PanelPanelInteractions.o		\
       TaylorMaster.o			\
       RWGGeometry.o			\
       Visualize.o			\
       RWGMisc.o			\

##################################################
##################################################
##################################################
       
##################################################
##################################################
##################################################
DEBUG=-O3
CC=g++
CPPFLAGS+=$(DEBUG) -I../SGJC -I../libTriInt
CPPFLAGS+=-I../libMatProp -I../libhmat -I../libSpherical
CPPFLAGS+=-I../libhrutil
CPPFLAGS+=$(GOOGLE_INCLUDES)

libscuff.a:	$(HEADERS) $(OBJS)
		ar r $@ $^

install:	libscuff.a
		(ln -s $(PWD)/libscuff.h $(HOME)/include)
		(ln -s $(PWD)/libscuff.a $(HOME)/lib)

clean:		
		rm -f *.o *.a
