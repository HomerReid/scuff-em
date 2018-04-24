ACLOCAL_AMFLAGS = -I m4

AM_CFLAGS = -O3
AM_CXXFLAGS = -O3
SUBDIRS = libs applications examples tests m4
EXTRA_DIST = NEWS.md LICENSE COPYRIGHT m4 scuff-em-pkgconfig.in

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = scuff-em.pc

scuff-em.pc: 	scuff-em-pkgconfig
		cp -f $(top_builddir)/scuff-em-pkgconfig $@

DISTCLEANFILES = scuff-em.pc
