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
STFT_OBJS = scuff-test-FIPPITable.o
STQDFIPPI_OBJS = scuff-test-QDFIPPIs.o GetQDFIPPIs_BruteForce.o
STQIFIPPI_OBJS = scuff-test-QIFIPPIs.o GetQIFIPPIs_BruteForce.o
STPPI_OBJS = scuff-test-PPIs.o GetPPIs_BruteForce.o

scuff-test-FIPPITable:	$(STFT_OBJS) libscuff.a
			$(CXX) $(LDFLAGS) -o $@ $(STFT_OBJS) $(LIBS)


scuff-test-PPIs:	$(STPPI_OBJS) libscuff.a
			$(CXX) $(LDFLAGS) -o $@ $(STPPI_OBJS) $(LIBS)

scuff-test-QDFIPPIs:	$(STQDFIPPI_OBJS) libscuff.a
			$(CXX) $(LDFLAGS) -o $@ $(STQDFIPPI_OBJS) $(LIBS)

scuff-test-QIFIPPIs:	$(STQIFIPPI_OBJS) libscuff.a
			$(CXX) $(LDFLAGS) -o $@ $(STQIFIPPI_OBJS) $(LIBS)

libscuff.a:	
		(pushd ../src/libs/libscuff; make DEBUG=$(DEBUG) libscuff.a; popd)

clean:	
		rm -f *.o *.a
		(pushd ../src/libs/libscuff; make clean; popd)
