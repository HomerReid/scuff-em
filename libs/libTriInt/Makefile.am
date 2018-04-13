noinst_LTLIBRARIES = libTriInt.la
pkginclude_HEADERS = libTriInt.h BZIntegration.h
libTriInt_la_SOURCES =  	\
 libTriInt.cc			\
 CCRules.cc			\
 LebedevRules.cc		\
 CCCubature.cc			\
 TCR.cc 			\
 BZIntegration.cc		\
 LatticeSum.cc   		\
 IntegrateCliffFunction.cc

AM_CPPFLAGS = -I$(top_srcdir)/libs/libhrutil \
              -I$(top_srcdir)/libs/libSGJC   \
              -I$(top_srcdir)/libs/libhmat 
