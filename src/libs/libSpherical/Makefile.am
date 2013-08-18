noinst_LTLIBRARIES = libSpherical.la
pkginclude_HEADERS = libSpherical.h AmosBessel.h
libSpherical_la_SOURCES = \
 libSpherical.h		\
 AmosBessel.h		\
 libSpherical.cc 	\
 AmosBessel.cc   	\
 TranslationMatrices.cc	\
 zbsubs.c		\
 machcon.c		\
 drc3jm.c 		

AM_CPPFLAGS = -I$(top_srcdir)/src/libs/libhrutil \
              -I$(top_srcdir)/src/libs/libhmat
