# Copyright(c)'1994-2009 by The Givaro group
# This file is part of Givaro.
# Givaro is governed by the CeCILL-B license under French law
# and abiding by the rules of distribution of free software.
# see the COPYRIGHT file for more details.

AM_CPPFLAGS=-I$(top_srcdir)  ${GMP_VERSION}

AM_CXXFLAGS = @DEFAULT_CFLAGS@
AM_CPPFLAGS+= -I$(top_srcdir)/src/kernel/integer -I$(top_srcdir)/src/kernel -I$(top_srcdir)/src/kernel/system -I$(top_srcdir)/src/kernel/memory -I$(top_srcdir)/src/kernel/bstruct

pkginclude_HEADERS=     \
	StaticElement.h		\
	gf2.h               \
	gfq.h		    	\
	gfqext.h		    \
	gfqkronecker.h		\
	extension.h			\
	givprimes16.h		\
	givrns.h		    \
	givrnsfixed.h		\
	chineseremainder.h	\
	givtablelimits.h	\
	givrnsfixed.inl		\
	givrnscstor.inl		\
	givrnsconvert.inl	\
	gf2.inl             \
	gfq.inl

noinst_LTLIBRARIES=libgivfield.la

libgivfield_la_SOURCES=   \
	givprimes16.C		\
	givtablelimits.C

EXTRA_DIST=field.doxy
