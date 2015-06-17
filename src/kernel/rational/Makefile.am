# Copyright(c)'1994-2009 by The Givaro group
# This file is part of Givaro.
# Givaro is governed by the CeCILL-B license under French law
# and abiding by the rules of distribution of free software.
# see the COPYRIGHT file for more details.

AM_CPPFLAGS=-I$(top_srcdir)  ${GMP_VERSION}

AM_CXXFLAGS = @DEFAULT_CFLAGS@
AM_CPPFLAGS+= -I$(top_srcdir)/src/kernel/integer -I$(top_srcdir)/src/kernel -I$(top_srcdir)/src/kernel/system -I$(top_srcdir)/src/kernel/memory -I$(top_srcdir)/src/kernel/ring

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
