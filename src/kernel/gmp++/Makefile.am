# Copyright(c)'1994-2009 by The Givaro group
# This file is part of Givaro.
# Givaro is governed by the CeCILL-B license under French law
# and abiding by the rules of distribution of free software.
# see the COPYRIGHT file for more details.

AM_CPPFLAGS=-I$(top_builddir) ${GMP_VERSION}

pkgincludesubdir=$(includedir)/gmp++

AM_CXXFLAGS=$(GIVARO_CXXFLAGS)
AM_CPPFLAGS+= -I$(top_builddir)/src/kernel/memory -I$(top_builddir)/src/kernel/system -I$(top_builddir)/src/kernel/recint


noinst_LTLIBRARIES=libgmppp.la

pkgincludesub_HEADERS=\
	gmp++.h       \
	gmp++_int.h   \
	gmp++_int_rand.inl

libgmppp_la_SOURCES= gmp++_int_lib.C

COMMONFILES= \
	gmp++_int_div.C		\
	gmp++_int_mod.C		\
	gmp++_int_sub.C		\
	gmp++_int_compare.C	\
	gmp++_int_gcd.C		\
	gmp++_int_io.C		\
	gmp++_int_mul.C		\
	gmp++_int_add.C		\
	gmp++_int_cstor.C	\
	gmp++_int_misc.C	\
	gmp++_int_pow.C


if GIVARO_INLINE_ALL

pkgincludesub_HEADERS += \
	gmp++_int.C             \
	$(COMMONFILES)


else

libgmppp_la_SOURCES+=            \
	$(COMMONFILES)

endif


EXTRA_DIST= gmp++.doxy
