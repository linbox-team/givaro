# Copyright(c)'1994-2009 by The Givaro group
# This file is part of Givaro.
# Givaro is governed by the CeCILL-B license under French law
# and abiding by the rules of distribution of free software.
# see the COPYRIGHT file for more details.

AM_CPPFLAGS=-I$(top_builddir)  ${GMP_VERSION}

AM_CXXFLAGS = $(GIVARO_CXXFLAGS)
AM_CPPFLAGS+= -I$(top_builddir)/src/kernel/integer -I$(top_builddir)/src/kernel -I$(top_builddir)/src/kernel/system -I$(top_builddir)/src/kernel/memory -I$(top_builddir)/src/kernel/bstruct

pkginclude_HEADERS=     \
	modular.h		    \
	modular-general.h	\
        modular-implem.h	\
	modular-integral.h	\
	modular-floating.h	\
	modular-integer.h	\
	modular-inttype.h	\
	modular-log16.h		\
	modular-mulprecomp.inl  \
	modular-ruint.h		\
	modular-balanced.h 			\
	modular-balanced-int32.h 	\
	modular-balanced-int64.h 	\
	modular-balanced-float.h 	\
	modular-balanced-double.h 	\
	unparametric-operations.h	\
	zring.h				\
	montgomery.h	    \
	montgomery-int32.h	\
	montgomery-ruint.h	\
	modular-general.inl	\
        modular-integral.inl    \
	modular-floating.inl	\
	modular-integer.inl	\
	modular-inttype.inl	\
	modular-log16.inl	\
	modular-ruint.inl   \
	modular-balanced-int32.inl  \
	modular-balanced-int64.inl 	\
	modular-balanced-float.inl 	\
	modular-balanced-double.inl \
	montgomery-int32.inl\
	montgomery-ruint.inl\
	ring-interface.h        \
	modular-extended.h	\
	modular-extended.inl

EXTRA_DIST=ring.doxy
