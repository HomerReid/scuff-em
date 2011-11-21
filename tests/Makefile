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

CPPFLAGS+=-I../src/libscuff
LDFLAGS+=-L../src/libscuff

##################################################
##################################################
##################################################
STPPI_OBJS = scuff-test-PPIs.o GetPPIs_BruteForce.o

scuff-test-PPIs:	$(STPPI_OBJS) libscuff.a
			$(CXX) $(LDFLAGS) -o $@ $(STPPI_OBJS) $(LIBS)

libscuff.a:	
		(pushd ../src/libscuff; make DEBUG=$(DEBUG) libscuff.a; popd)

clean:	
		rm -f *.o *.a
		(pushd ../src/libscuff; make clean; popd)
