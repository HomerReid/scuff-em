######################################################################
# include system-dependent paths etc.
#######################################################################
include $(HOME)/work/Makefile.Rules

##################################################
##################################################
##################################################
HEADERS = libscuff.h			\
          libscuffInternals.h

OBJS = AssembleBEMMatrix.o		\
       AssembleDMDVMatrix.o		\
       AssembleRHSVector.o		\
       CalcGC.o				\
       GetDipoleMoments.o		\
       GetSphericalMoments.o		\
       GetFields.o			\
       GetMatrixElement.o		\
       ExpandCurrentDistribution.o	\
       InitEdgeList.o			\
       RWGObject.o			\
       LFunctions.o			\
       ReadComsolFile.o			\
       ReadGMSHFile.o			\
       StaticPPI.o			\
       StaticPPITaylor.o		\
       TaylorMaster.o			\
       RWGGeometry.o			\
       Visualize.o			\
       RWGMisc.o			\
       AssembleRHSVectorArray.o		\
       GetFieldsArray.o        		\
       GetStressTensorArray.o

##################################################
##################################################
##################################################
       
##################################################
##################################################
##################################################
DEBUG=-O3
CC=g++
CPPFLAGS+=$(DEBUG) -I../SGJV -I../libTriInt
CPPFLAGS+=-I../libMatProp -I../libhmat -I../libSpherical
CPPFLAGS+=-I../libhrutil
CPPFLAGS+=$(GOOGLE_INCLUDES)

libRWG.a:	$(HEADERS) $(OBJS)
		ar r $@ $^

install:	libRWG.a
		(ln -s $(PWD)/libRWG.h $(HOME)/include)
		(ln -s $(PWD)/libRWG.a $(HOME)/lib)

clean:		
		rm -f *.o *.a
