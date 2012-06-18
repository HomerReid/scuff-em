#
# this is a simple, handwritten, bare-bones makefile that 
# indicates how to compile and link a C++ program with 
# libscuff.
#
# homer reid 6/2012
#

# 
# root directory of scuff-EM installation 
#

SCUFFEM = $(HOME)/scuff-em-installation
#SCUFFEM = /usr/local

CPPFLAGS += -I$(SCUFFEM)/include
LDFLAGS += -L$(SCUFFEM)/lib

#
# you may need to change this line depending on the 
# names of the HDF5 and lapack/blas library binaries 
# on your system
#

LIBS = -lscuff -lhdf5_hl -llapack -lblas

#
# hopefully this should be all that is needed
#
PlateWithHole:	PlateWithHole.o
		$(CXX) -o PlateWithHole PlateWithHole.o $(LIBS)
