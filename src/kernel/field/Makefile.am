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
	givgfq.h		    \
	modular.h		    \
	modular-int16.h		\
	modular-int32.h		\
	modular-int64.h		\
	modular-uint32.h	\
	modular-integer.h	\
	modular-inttype.h	\
	modular-log16.h		\
	montgomery.h	    \
	montgomery-int32.h	\
	givprimes16.h		\
	givrns.h		    \
	givrnsfixed.h		\
	givcra.h		    \
	giv_randiter.h		\
	givextension.h		\
	givgfqext.h		    \
	givtablelimits.h	\
	givgfqkron.h		\
	modular-int16.inl	\
	modular-int32.inl	\
	modular-int64.inl	\
	modular-uint32.inl	\
	modular-integer.inl	\
	modular-inttype.inl	\
	modular-log16.inl	\
	montgomery-int32.inl\
	givrnsfixed.inl		\
	givrnscstor.inl		\
	givrnsconvert.inl	\
	givgfq.inl

noinst_LTLIBRARIES=libgivfield.la

libgivfield_la_SOURCES=   \
	modular-int16.C		\
	modular-int32.C		\
	modular-int64.C		\
	modular-uint32.C	\
	modular-log16.C		\
	montgomery-int32.C	\
	givprimes16.C		\
	givtablelimits.C

EXTRA_DIST=field.doxy
