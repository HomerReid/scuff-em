# Copyright (C) 1999, 2002, 2003, 2004, 2005, 2006, 2007 Free Software
# Foundation, Inc.
#
# This file is part of GNU libmatheval
#
# GNU libmatheval is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License as
# published by the Free Software Foundation; either version 2, or (at
# your option) any later version.
#
# GNU libmatheval is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with program; see the file COPYING. If not, write to the Free
# Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
# 02110-1301, USA.

CC = $(CXX)
AM_CFLAGS = $(CXXFLAGS)

MAINTAINERCLEANFILES = Makefile.in parser.c parser.h scanner.c

AM_CPFLAGS = -I$(top_builddir) -I$(top_srcdir)

noinst_LTLIBRARIES = libcmatheval.la

libcmatheval_la_SOURCES = parser.y scanner.l matheval.c	\
node.c symbol_table.c xmalloc.c xmath.c

pkginclude_HEADERS = cmatheval.h
noinst_HEADERS = common.h node.h symbol_table.h xmalloc.h xmath.h

AM_YFLAGS = -d

noinst_PROGRAMS = tcmatheval
tcmatheval_SOURCES = tcmatheval.c
tcmatheval_LDADD = libcmatheval.la -lm
