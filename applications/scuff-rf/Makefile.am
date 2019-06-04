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

scuff_rf_LDADD = $(top_builddir)/libs/libscuff/libscuff.la

AM_CPPFLAGS = -DSCUFF \
              -I$(top_srcdir)/libs/libhrutil	\
              -I$(top_srcdir)/libs/libscuff		 \
              -I$(top_srcdir)/libs/libIncField   \
              -I$(top_srcdir)/libs/libMatProp    \
              -I$(top_srcdir)/libs/libMDInterp    \
              -I$(top_srcdir)/libs/libhmat       \
              -I$(top_srcdir)/libs/libSGJC       \
              -I$(top_srcdir)/libs/libSubstrate  \
              -I$(top_srcdir)/libs/libTriInt     \
              -DVERSIONSTRING=\"$(shell git describe always --dirty)\"
