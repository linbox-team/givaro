# Copyright(c)'1994-2009 by The Givaro group
# This file is part of Givaro.
# Givaro is governed by the CeCILL-B license under French law
# and abiding by the rules of distribution of free software.
# see the COPYRIGHT file for more details.

INCLUDES=-I$(top_srcdir)

AM_CXXFLAGS=@DEFAULT_CFLAGS@
AM_CPPFLAGS= -I$(top_srcdir)/src/kernel/system

pkginclude_HEADERS= \
	givaromm.h     	\
	givpointer.h   	\
	givref_count.h
noinst_LTLIBRARIES=libgivmemory.la

libgivmemory_la_SOURCES=givaromm.C
