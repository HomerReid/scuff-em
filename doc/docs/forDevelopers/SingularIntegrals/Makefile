# top of SCUFF-EM src tree
SCUFF_SRC=/home/homer/work/scuff-em

# SCUFF-EM installation prefix (set with --prefix when running configure)
SCUFF_INSTALL=/home/homer/work/scuff-em-installation

# HDF5 / Lapack libraries
HDF5_LIBS=-lhdf5_hl -lhdf5
MATH_LIBS=-lopenblas -lgomp -lgfortran -lpthread

##################################################
# shouldn't need to modify what follows
##################################################
CPPFLAGS+=-I$(SCUFF_SRC)/src/libs/libscuff 
CPPFLAGS+=-I$(SCUFF_INSTALL)/include/scuff-em

LIBDIR=$(SCUFF_INSTALL)/lib
LDFLAGS+=-L$(LIBDIR) -Wl,-rpath,$(LIBDIR)

SCUFF_LIBS=-lscuff
LIBS=$(SCUFF_LIBS) $(HDF5_LIBS) $(MATH_LIBS)

TestTaylorDuffy:		TestTaylorDuffy.o
				$(CXX) $(LDFLAGS) -o $@ $^ $(LIBS) 
