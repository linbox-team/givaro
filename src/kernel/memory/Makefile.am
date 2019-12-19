# Copyright(c)'1994-2009 by The Givaro group
# This file is part of Givaro.
# Givaro is governed by the CeCILL-B license under French law
# and abiding by the rules of distribution of free software.
# see the COPYRIGHT file for more details.

AM_CPPFLAGS=-I$(top_builddir)

AM_CXXFLAGS=$(GIVARO_CXXFLAGS)
AM_CPPFLAGS+= -I$(top_builddir)/src/kernel/system

pkginclude_HEADERS= \
	givaromm.h     	\
	givmemorer.h     	\
	givpointer.h   	\
	givref_count.h
noinst_LTLIBRARIES=libgivmemory.la

libgivmemory_la_SOURCES=givaromm.C

EXTRA_DIST=memory.doxy
