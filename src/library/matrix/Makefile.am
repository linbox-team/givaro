# Copyright(c)'1994-2009 by The Givaro group
# This file is part of Givaro.
# Givaro is governed by the CeCILL-B license under French law
# and abiding by the rules of distribution of free software.
# see the COPYRIGHT file for more details.

AM_CPPFLAGS=-I$(top_builddir)

AM_CXXFLAGS = $(GIVARO_CXXFLAGS)
AM_CPPFLAGS+= -I$(top_builddir)/src/kernel/bstruct -I$(top_builddir)/src/kernel/system -I$(top_builddir)/src/kernel/memory -I$(top_builddir)/src/library/vector

pkginclude_HEADERS=         \
	givmatdense.h		    \
	givmatrix.h		        \
	givmatsparse.h		    \
	givmatstoragedense.h	\
	givmatsparseops.inl	    \
	givmatstoragesparse.h	\
	givmatdenseops.inl	    \
	givmatstorage.h

#noinst_LTLIBRARIES=libgivmatrix.la
#libgivmatrix_la_SOURCES=givmatrix.C
