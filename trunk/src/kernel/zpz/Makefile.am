# Copyright(c)'1994-2009 by The Givaro group
# This file is part of Givaro.
# Givaro is governed by the CeCILL-B license under French law
# and abiding by the rules of distribution of free software.
# see the COPYRIGHT file for more details.

INCLUDES=-I$(top_srcdir)  ${GMP_VERSION}

AM_CXXFLAGS = @DEFAULT_CFLAGS@
AM_CPPFLAGS= -I$(top_srcdir)/src/kernel/integer -I$(top_srcdir)/src/kernel -I$(top_srcdir)/src/kernel/system -I$(top_srcdir)/src/kernel/memory -I$(top_srcdir)/src/kernel/bstruct

pkginclude_HEADERS=\
	StaticElement.h		\
	givzpz.h		\
	givzpztypes.h		\
	givgfq.h		\
	givzpz16std.h		\
	givzpz32std.h		\
	givzpz32uns.h		\
	givzpzInt.h		\
	givzpz64std.h		\
	givprimes16.h		\
	givrns.h		\
	givcra.h		\
	givzpz16table1.h	\
	givmontg32.h		\
	giv_randiter.h		\
	givextension.h		\
	givgfqext.h		\
	givtablelimits.h	\
	givgfqkron.h		\
	givzpz16std.inl		\
	givzpz32std.inl		\
	givzpz32uns.inl		\
	givzpzInt.inl		\
	givzpz64std.inl		\
	givrnscstor.inl		\
	givrnsconvert.inl	\
	givzpz16table1.inl	\
	givmontg32.inl		\
	givgfq.inl

noinst_LTLIBRARIES=libgivzpz.la

libgivzpz_la_SOURCES=\
	givzpz16std.C		\
	givzpz32std.C		\
	givzpz32uns.C		\
	givmontg32.C		\
	givzpzInt.C		\
	givzpz64std.C		\
	givprimes16.C		\
	givtablelimits.C	\
	givzpz16table1.C
