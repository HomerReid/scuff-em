#
# NOTE: the source files 'SSSolver.cc', 'GetPPI.cc',
# 'GetPhiE.cc,' and 'Substrate.cc' together constitute the 
# the full implementation of the SSSolver class. The files 
# 'scuff-static.cc' and  'OutputModules.cc' are just one 
# particular instance of a driver program that uses this class. 
# In principle the former files should be compiled into the form
# of a library, which is then linked by the scuff-static 
# executable but could be linked by other executables too.
# But here for convenience we just throw them all into one 
# big executable.
# 

bin_PROGRAMS = scuff-static

scuff_static_SOURCES = 		\
 GetPhiE.cc			\
 GetPPI.cc			\
 SSSolver.cc			\
 SSSolver.h			\
 OutputModules.cc		\
 StaticSubstrate.h     		\
 StaticSubstrate.cc    		\
 PanelPanelCubature2.cc		\
 scuff-static.cc

scuff_static_LDADD = $(top_builddir)/src/libs/libscuff/libscuff.la

AM_CPPFLAGS = -I$(top_srcdir)/src/libs/libscuff      \
              -I$(top_srcdir)/src/libs/libIncField   \
              -I$(top_srcdir)/src/libs/libMatProp    \
              -I$(top_srcdir)/src/libs/libMatProp/cmatheval \
              -I$(top_srcdir)/src/libs/libMDInterp   \
              -I$(top_srcdir)/src/libs/libhmat       \
              -I$(top_srcdir)/src/libs/libSGJC       \
              -I$(top_srcdir)/src/libs/libSpherical  \
              -I$(top_srcdir)/src/libs/libTriInt     \
              -I$(top_srcdir)/src/libs/libhrutil
