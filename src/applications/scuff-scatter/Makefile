SHELL=/bin/bash

##################################################
##################################################
##################################################
include ../../Makefile.Rules

##################################################
##################################################
##################################################
HR_LIBS = -lscuff -lPolyFit -lMatProp -lMDInterp \
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

CPPFLAGS+=-I../src/libs/libscuff
LDFLAGS+=-L../src//libs/libscuff

##################################################
##################################################
##################################################
SS_OBJS = scuff-scatter.o OutputModules.o

scuff-scatter:	$(SS_OBJS) libscuff.a
		$(CXX) $(LDFLAGS) -o $@ $(SS_OBJS) $(LIBS)

libscuff.a:	
		(pushd ../src/libs/libscuff; make DEBUG=$(DEBUG) libscuff.a; popd)

clean:	
		rm -f *.o *.a
		(pushd ../src/libs/libscuff; make clean; popd)
