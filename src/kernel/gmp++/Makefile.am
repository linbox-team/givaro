# Copyright(c)'1994-2009 by The Givaro group
# This file is part of Givaro.
# Givaro is governed by the CeCILL-B license under French law
# and abiding by the rules of distribution of free software.
# see the COPYRIGHT file for more details.

INCLUDES=-I$(top_srcdir) ${GMP_VERSION}

pkgincludesubdir=$(includedir)/gmp++

AM_CXXFLAGS=@DEFAULT_CFLAGS@
AM_CPPFLAGS= -I$(top_srcdir)/src/kernel/memory -I$(top_srcdir)/src/kernel/system

pkgincludesub_HEADERS=\
	gmp++.h		\
	gmp++_int.h	\
	gmp++_int.inl	\
	gmp++.C		\
	gmp++_int.C

#lib_LTLIBRARIES=libgmpxx.la
noinst_LTLIBRARIES=libgmppp.la

#libgmpxx_la_SOURCES=
libgmppp_la_SOURCES=\
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

#EXTRA_DIST= gmp++.doxy
