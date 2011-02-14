# Copyright(c)'1994-2009 by The Givaro group
# This file is part of Givaro.
# Givaro is governed by the CeCILL-B license under French law
# and abiding by the rules of distribution of free software.
# see the COPYRIGHT file for more details.

INCLUDES=-I$(top_srcdir)

AM_CXXFLAGS = @DEFAULT_CFLAGS@
AM_CPPFLAGS= -I$(top_srcdir)/src/kernel/bstruct -I$(top_srcdir)/src/kernel/system -I$(top_srcdir)/src/kernel/memory -I$(top_srcdir)/src/library/vector

pkginclude_HEADERS=\
	givmatdense.h		\
	givmatrix.h		\
	givmatsparse.h		\
	givmatstoragedense.h	\
	givmatsparseops.inl	\
	givmatstoragesparse.h	\
	givmatdenseops.inl	\
	givmatstorage.h

noinst_LTLIBRARIES=libgivmatrix.la

libgivmatrix_la_SOURCES=givmatrix.C
