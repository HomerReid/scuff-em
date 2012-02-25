bin_PROGRAMS = scuff-scatter

scuff_scatter_SOURCES = 	\
 scuff-scatter.cc 		\
 OutputModules.cc 		\
 scuff-heat.h

LDADD = $(top_srcdir)/src/libs/libscuff/libscuff.a	\
        $(top_srcdir)/src/libs/libIncField/libIncField.a  \
        $(top_srcdir)/src/libs/libMatProp/libMatProp.a	\
        $(top_srcdir)/src/libs/libhmat/libhmat.a      	\
        $(top_srcdir)/src/libs/libSGJC/libSGJC.a      	\
        $(top_srcdir)/src/libs/libTriInt/libTriInt.a    \
        $(top_srcdir)/src/libs/libTriInt/libMDInterp.a  \
        $(top_srcdir)/src/libs/libMDInterp/libhrutil.a

AM_CPPFLAGS = -I$(top_srcdir)/src/libs/libscuff      \
              -I$(top_srcdir)/src/libs/libIncField   \
              -I$(top_srcdir)/src/libs/libMatProp    \
              -I$(top_srcdir)/src/libs/libhmat       \
              -I$(top_srcdir)/src/libs/libSGJC       \
              -I$(top_srcdir)/src/libs/libTriInt     \
              -I$(top_srcdir)/src/libs/libhrutil
