SHELL=/bin/bash

##################################################
##################################################
##################################################
include ../../../../Makefile.Rules

##################################################
##################################################
##################################################
#HR_LIBS = -lscuff -lPolyFit -lMatProp -lMDInterp \
#          -lhmat -lSGJC -lTriInt -lhrutil

HR_LIBS = -lRWG -lIncFieldPrime -lMatProp -lMDInterp \
          -lhmat -lSGJC -lTriInt -lhrutil

OTHER_LIBS = $(GSL_LIBS) $(MUPARSER_LIBS) $(RDL_LIBS) $(HDF5_LIBS) $(LB_LIBS)

LIBS=$(HR_LIBS) $(OTHER_LIBS)

##################################################
##################################################
##################################################
DEBUG=-O2
CXX=g++
CFLAGS+=$(DEBUG)
CPPFLAGS+=$(DEBUG)
LDFLAGS+=$(DEBUG)

CPPFLAGS+=-I../../libs/libscuff
CPPFLAGS+=-I../../libs/libIncField

LDFLAGS+=-L.
LDFLAGS+=-L../../libs/libscuff
LDFLAGS+=-L../../libs/libIncField

##################################################
##################################################
##################################################
SS_OBJS = scuff-scatter.o OutputModules.o

scuff-scatter:	$(SS_OBJS) 
		$(CXX) $(LDFLAGS) -o $@ $(SS_OBJS) $(LIBS)

#libscuff.a:	
#		(pushd ../src/libs/libscuff; make DEBUG=$(DEBUG) libscuff.a; popd)

clean:	
		rm -f *.o *.a
#		(pushd ../src/libs/libscuff; make clean; popd)
