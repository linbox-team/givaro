# Copyright(c)'1994-2009 by The Givaro group
# This file is part of Givaro.
# Givaro is governed by the CeCILL-B license under French law
# and abiding by the rules of distribution of free software.
# see the COPYRIGHT file for more details.

AM_CPPFLAGS=-I$(top_srcdir) ${GMP_VERSION}

pkgincludesubdir=$(includedir)/gmp++

AM_CXXFLAGS=@DEFAULT_CFLAGS@
AM_CPPFLAGS+= -I$(top_srcdir)/src/kernel/memory -I$(top_srcdir)/src/kernel/system


noinst_LTLIBRARIES=libgmppp.la

pkgincludesub_HEADERS=\
	gmp++.h       \
	gmp++_int.h   \
	gmp++_rat.h   \
	gmp++_rat.inl

libgmppp_la_SOURCES= gmp++_int_lib.C

if GIVARO_INLINE_ALL

pkgincludesub_HEADERS += \
	gmp++_int.C             \
	gmp++_int_rand.inl	\
	gmp++_rat.C             \
	gmp++_int_div.C		\
	gmp++_int_mod.C		\
	gmp++_int_sub.C		\
	gmp++_int_compare.C	\
	gmp++_rat_compare.C	\
	gmp++_int_gcd.C		\
	gmp++_int_io.C		\
	gmp++_rat_io.C		\
	gmp++_int_mul.C		\
	gmp++_int_add.C		\
	gmp++_rat_add.C		\
	gmp++_int_cstor.C	\
	gmp++_rat_cstor.C	\
	gmp++_int_misc.C	\
	gmp++_rat_misc.C	\
	gmp++_int_pow.C


else

libgmppp_la_SOURCES+=            \
	gmp++_int_div.C		\
	gmp++_int_mod.C		\
	gmp++_int_sub.C		\
	gmp++_int_compare.C	\
	gmp++_rat_compare.C	\
	gmp++_int_gcd.C		\
	gmp++_int_io.C		\
	gmp++_rat_io.C		\
	gmp++_int_mul.C		\
	gmp++_int_add.C		\
	gmp++_rat_add.C		\
	gmp++_int_cstor.C	\
	gmp++_rat_cstor.C	\
	gmp++_int_misc.C	\
	gmp++_rat_misc.C	\
	gmp++_int_pow.C      \
	gmp++_int_lib.C

endif


EXTRA_DIST= gmp++.doxy
