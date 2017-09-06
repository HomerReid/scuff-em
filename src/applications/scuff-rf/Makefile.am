bin_PROGRAMS = scuff-rf

scuff_rf_SOURCES =		\
 scuff-rf.cc			\
 EdgePanelInteractions.cc	\
 GetPanelPotentials.cc		\
 GetPortVoltages.cc		\
 OutputModules.cc		\
 RWGPorts.cc			\
 RWGPorts.h 	 		\
 ZSConvert.cc

scuff_rf_LDADD = $(top_builddir)/src/libs/libscuff/libscuff.la

AM_CPPFLAGS = -DSCUFF \
              -I$(top_srcdir)/src/libs/libhrutil	\
              -I$(top_srcdir)/src/libs/libscuff		 \
              -I$(top_srcdir)/src/libs/libIncField   \
              -I$(top_srcdir)/src/libs/libMatProp    \
              -I$(top_srcdir)/src/libs/libMDInterp    \
              -I$(top_srcdir)/src/libs/libhmat       \
              -I$(top_srcdir)/src/libs/libSGJC       \
              -I$(top_srcdir)/src/libs/libTriInt
