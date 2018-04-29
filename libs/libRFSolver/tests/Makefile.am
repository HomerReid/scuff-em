AM_CPPFLAGS = -I$(top_srcdir)/libs/libhrutil    \
              -I$(top_srcdir)/libs/libhmat      \
              -I$(top_srcdir)/libs/libMatProp   \
              -I$(top_srcdir)/libs/libMDInterp  \
              -I$(top_srcdir)/libs/libSGJC      \
              -I$(top_srcdir)/libs/libSpherical \
              -I$(top_srcdir)/libs/libIncField  \
              -I$(top_srcdir)/libs/libSubstrate \
              -I$(top_srcdir)/libs/libRFSolver  \
              -I$(top_srcdir)/libs/libscuff     \
              -I$(top_builddir) # for config.h

TESTLIBS = $(top_builddir)/libs/libscuff/libscuff.la

check_PROGRAMS = t1BFMOIFields tEquivalentPairs

t1BFMOIFields_SOURCES = t1BFMOIFields.cc
t1BFMOIFields_LDADD   = $(TESTLIBS)

tEquivalentPairs_SOURCES = tEquivalentPairs.cc
tEquivalentPairs_LDADD   = $(TESTLIBS)

 
#TESTS = t1BFMOIFields
