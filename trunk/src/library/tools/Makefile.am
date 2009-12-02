# Copyright(c)'1994-2009 by The Givaro group
# This file is part of Givaro.
# Givaro is governed by the CeCILL-B license under French law
# and abiding by the rules of distribution of free software. 
# see the COPYRIGHT file for more details.

INCLUDES=-I$(top_srcdir)

AM_CPPFLAGS= -I$(prefix)/include

pkginclude_HEADERS=\
	givcategory.h		\
	givops.h		\
	giviterator.h		\
	givarithmetics.h	\
	givfractiondomain.h	\
	givquotientdomain.h	\
	givtruncdomain.h	\
	givtruncdomain.inl	\
	givfixedtrunc.h		\
	givhighorder.h

noinst_LTLIBRARIES=libgivtools.la

libgivtools_la_SOURCES=givops.C
