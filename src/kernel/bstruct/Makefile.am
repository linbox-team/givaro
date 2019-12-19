# Copyright(c)'1994-2009 by The Givaro group
# This file is part of Givaro.
# Givaro is governed by the CeCILL-B license under French law
# and abiding by the rules of distribution of free software.
# see the COPYRIGHT file for more details.

AM_CPPFLAGS=-I$(top_builddir)

AM_CXXFLAGS = $(GIVARO_CXXFLAGS)
AM_CPPFLAGS+= -I$(top_builddir)/src/kernel/system -I$(top_builddir)/src/kernel/memory

pkginclude_HEADERS=     \
	givarray0.h		    \
	givarrayallocator.h	\
	givarrayfixed.h		\
	givbits.h		    \
	givelem.h		    \
	givhashtable.h		\
	givlist0.h		    \
	givstack.h		    \
	givarray0.inl		\
	givbits.inl		    \
	givhashtable.inl	\
	givlist0.inl		\
	givstack.inl

EXTRA_DIST = bstruct.doxy

noinst_LTLIBRARIES=libgivbstruct.la

libgivbstruct_la_SOURCES=givbits.C
