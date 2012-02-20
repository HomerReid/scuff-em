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

CPPFLAGS+=-DHAVE_MUPARSER

##################################################
##################################################
##################################################
TBI_OBJS = testB0Integrals.o

testB0Integrals:	$(TBI_OBJS)
			$(CXX) $(LDFLAGS) -o $@ $(TBI_OBJS) $(LIBS)

clean:	
		rm -f *.o *.a
