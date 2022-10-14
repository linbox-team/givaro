# Copyright(c)'1994-2009 by The Givaro group
# This file is part of Givaro.
# Givaro is governed by the CeCILL-B license under French law
# and abiding by the rules of distribution of free software.
# see the COPYRIGHT file for more details.

AM_CPPFLAGS=-I$(top_builddir)

AM_CXXFLAGS = $(GIVARO_CXXFLAGS)
AM_CPPFLAGS+= -I$(top_builddir)/src/kernel/bstruct -I$(top_builddir)/src/kernel/system -I$(top_builddir)/src/kernel/memory

pkginclude_HEADERS=        \
	givindeter.h		   \
	givpoly1.h			   \
	givdegree.h			   \
	givinterp.h			   \
	givinterpgeom.h		   \
	givinterpgeom-multip.h \
	givpoly1dense.h		   \
	givpoly1factor.h	   \
	givpoly1padic.h		   \
	givpoly1denseops.inl   \
	givpoly1addsub.inl	   \
	givpoly1muldiv.inl	   \
	givpoly1axpy.inl	   \
	givpoly1gcd.inl		   \
	givpoly1kara.inl	   \
	givpoly1midmul.inl	   \
	givpoly1cstor.inl	   \
	givpoly1io.inl		   \
	givpoly1misc.inl	   \
	givpoly1sqrfree.inl	   \
	givpoly1cyclo.inl	   \
	givpoly1factor.inl	   \
	givpoly1proot.inl	   \
	givpoly1ratrecon.inl   \
	givpoly1crt.h		   \
	givpoly1crtcstor.inl   \
	givpoly1crtconvert.inl

noinst_LTLIBRARIES=libgivpoly1.la

libgivpoly1_la_SOURCES=givindeter.C givdegree.C
