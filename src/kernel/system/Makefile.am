# Copyright(c)'1994-2009 by The Givaro group
# This file is part of Givaro.
# Givaro is governed by the CeCILL-B license under French law
# and abiding by the rules of distribution of free software.
# see the COPYRIGHT file for more details.

INCLUDES=-I$(top_srcdir)

AM_CXXFLAGS=@DEFAULT_CFLAGS@
AM_CPPFLAGS= -I$(top_srcdir)/src/kernel/memory

pkginclude_HEADERS= \
	givbasictype.h	\
	givconfig.h	    \
	giverror.h	    \
	givgenarith.h	\
	givinit.h	    \
	givmodule.h	    \
	givperf.h 	    \
	givpower.h	    \
	givrandom.h	    \
	givtimer.h

noinst_LTLIBRARIES=libgivsystem.la

libgivsystem_la_SOURCES= \
	givbasictype.C	     \
	giverror.C	         \
	givinit.C	         \
	givmodule.C	         \
	givtimer.C
