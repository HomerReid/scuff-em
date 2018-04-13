SUBDIRS = cmatheval

noinst_LTLIBRARIES = libMatProp.la
pkginclude_HEADERS = libMatProp.h
libMatProp_la_SOURCES = libMatProp.cc UserDefinedMaterials.cc
libMatProp_la_LIBADD = $(builddir)/cmatheval/libcmatheval.la

EXTRA_DIST = README matprop.dat.example

AM_CPPFLAGS = -I$(top_srcdir)/libs/libhrutil	\
              -I$(top_srcdir)/libs/libMatProp/cmatheval	\
              -I$(top_srcdir)/libs/libMDInterp	\
              -I$(top_srcdir)/libs/libhmat          \
	      -I$(top_builddir) # for config.h

#noinst_PROGRAMS = tlibMatProp
#tlibMatProp_SOURCES = tlibMatProp.cc
#tlibMatProp_LDADD = libMatProp.la ../libMDInterp/libMDInterp.la ../libhmat/libhmat.la ../libhrutil/libhrutil.la
