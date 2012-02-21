SHELL=/bin/bash

##################################################
##################################################
##################################################
include ../../../../Makefile.Rules

##################################################
##################################################
##################################################
HR_LIBS = -lscuff -lMatProp -lMDInterp \
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
CPPFLAGS+=-DSCUFF -DHAVE_MUPARSER -DHAVE_HDF5

LDFLAGS+=-L.
LDFLAGS+=-L../../libs/libscuff

##################################################
##################################################
##################################################
SH_OBJS = scuff-heat.o FrequencyIntegrand.o

scuff-heat:   	$(SH_OBJS) libscuff.a 
		$(CXX) $(LDFLAGS) -o $@ $(SH_OBJS) $(LIBS)

libscuff.a:	
		(pushd ../../libs/libscuff; make DEBUG=$(DEBUG) libscuff.a; popd)

clean:	
		rm -f *.o *.a
		(pushd ../../libs/libscuff; make clean; popd)
