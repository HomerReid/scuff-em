SWIG_SRC = scuff.i scuff-python.i
EXTRA_DIST = $(SWIG_SRC) 
BUILT_SOURCES = scuff-python.cpp scuff.py

_scuff_la_SOURCES = scuff-python.cpp
_scuff_la_LIBADD = ../libscuff/libscuff.la $(LIBS)
_scuff_la_LDFLAGS = -module -version-info @SHARED_VERSION_INFO@
_scuff_la_CPPFLAGS = $(PYTHON_INCLUDES) $(AM_CPPFLAGS)

if WITH_PYTHON
python_PYTHON = scuff.py
pyexec_LTLIBRARIES = _scuff.la
endif

AM_CPPFLAGS = -I$(top_srcdir)/libs/libhrutil       \
              -I$(top_srcdir)/libs/libMatProp      \
              -I$(top_srcdir)/libs/libhmat         \
              -I$(top_srcdir)/libs/libSGJC         \
              -I$(top_srcdir)/libs/libIncField     \
              -I$(top_srcdir)/libs/libMDInterp     \
              -I$(top_srcdir)/libs/libTriInt       \
              -I$(top_srcdir)/libs/libscuff        \
              -I$(top_srcdir)/libs/libscuffSolver  \
              -I$(top_srcdir)/libs/libStaticSolver \
              -I$(top_srcdir)/libs/libSubstrate    \
              -I$(top_builddir) # for config.h

AM_CFLAGS=-W
AM_CXXFLAGS=-W

##############################################################################
# Generate wrappers via SWIG:

if MAINTAINER_MODE

scuff-python.cpp scuff.py: $(SWIG_SRC) $(top_srcdir)/libs/libscuff/libscuff.h $(top_srcdir)/libs/libscuff/GTransformation.h $(top_srcdir)/libs/libIncField/libIncField.h $(top_srcdir)/libs/libhrutil/libhrutil.h $(top_srcdir)/libs/libMatProp/libMatProp.h $(top_srcdir)/libs/libhmat/libhmat.h $(top_srcdir)/libs/libStaticSolver/StaticSolver.h $(top_srcdir)/libs/libscuffSolver/scuffSolver.h
	swig -Wall $(AM_CPPFLAGS) -outdir $(builddir) -c++ -python -o $(builddir)/scuff-python.cpp $(srcdir)/scuff.i

endif
