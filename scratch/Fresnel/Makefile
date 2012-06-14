SCUFF_INSTALLATION=/home/homer/work/scuff-em-installation

CXX = g++

CFLAGS = -O2 
CPPFLAGS = -I$(SCUFF_INSTALLATION)/include/scuff-em
LDFLAGS = -L$(SCUFF_INSTALLATION)/lib

LIBS = -lscuff -lhdf5 -lhdf5_hl -llapack -lblas -lm -lpthread -lgfortran
 
Fresnel:	Fresnel.o 
		$(CXX) $(LDFLAGS) -o $@ $^ $(LIBS) 
