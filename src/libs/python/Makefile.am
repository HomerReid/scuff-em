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

AM_CPPFLAGS = -I$(top_srcdir)/src/libs/libhrutil    \
              -I$(top_srcdir)/src/libs/libMatProp   \
              -I$(top_srcdir)/src/libs/libhmat      \
              -I$(top_srcdir)/src/libs/libSGJC      \
              -I$(top_srcdir)/src/libs/libIncField  \
              -I$(top_srcdir)/src/libs/libMDInterp  \
              -I$(top_srcdir)/src/libs/libTriInt    \
              -I$(top_srcdir)/src/libs/libscuff     \
              -I$(top_builddir) # for config.h

##############################################################################
# Generate wrappers via SWIG:

if MAINTAINER_MODE

scuff-python.cpp scuff.py: $(SWIG_SRC) $(top_srcdir)/src/libs/libscuff/libscuff.h $(top_srcdir)/src/libs/libscuff/GTransformation.h $(top_srcdir)/src/libs/libIncField/libIncField.h $(top_srcdir)/src/libs/libhrutil/libhrutil.h $(top_srcdir)/src/libs/libMatProp/libMatProp.h $(top_srcdir)/src/libs/libhmat/libhmat.h
	swig -Wall $(AM_CPPFLAGS) -outdir $(builddir) -c++ -python -o $(builddir)/scuff-python.cpp $(srcdir)/scuff.i

endif

