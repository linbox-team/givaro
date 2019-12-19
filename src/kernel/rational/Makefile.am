# Copyright(c)'1994-2009 by The Givaro group
# This file is part of Givaro.
# Givaro is governed by the CeCILL-B license under French law
# and abiding by the rules of distribution of free software.
# see the COPYRIGHT file for more details.

AM_CPPFLAGS=-I$(top_builddir)  ${GMP_VERSION}

AM_CXXFLAGS = $(GIVARO_CXXFLAGS)
AM_CPPFLAGS+= -I$(top_builddir)/src/kernel/integer -I$(top_builddir)/src/kernel -I$(top_builddir)/src/kernel/system -I$(top_builddir)/src/kernel/memory -I$(top_builddir)/src/kernel/ring -I$(top_srcdir)/src/kernel

pkginclude_HEADERS= \
	givrational.h	\
	givrational.inl

noinst_LTLIBRARIES=libgivrational.la

libgivrational_la_SOURCES= \
	givratcstor.C	       \
	givratio.C	           \
	givratmisc.C	       \
	givratmuldiv.C	       \
	givrataddsub.C	       \
	givratcpy.C	           \
	givratcompare.C        \
	givratreconstruct.C

EXTRA_DIST=rational.doxy
