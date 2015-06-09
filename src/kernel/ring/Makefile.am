# Copyright(c)'1994-2009 by The Givaro group
# This file is part of Givaro.
# Givaro is governed by the CeCILL-B license under French law
# and abiding by the rules of distribution of free software.
# see the COPYRIGHT file for more details.

AM_CPPFLAGS=-I$(top_srcdir)  ${GMP_VERSION}

AM_CXXFLAGS = @DEFAULT_CFLAGS@
AM_CPPFLAGS+= -I$(top_srcdir)/src/kernel/integer -I$(top_srcdir)/src/kernel -I$(top_srcdir)/src/kernel/system -I$(top_srcdir)/src/kernel/memory -I$(top_srcdir)/src/kernel/bstruct

pkginclude_HEADERS=     \
	modular.h		    \
	modular-general.h	\
	modular-defines.h	\
	modular-int8.h		\
	modular-int16.h		\
	modular-int32.h		\
	modular-int64.h		\
	modular-uint8.h		\
	modular-uint16.h	\
	modular-uint32.h	\
	modular-uint64.h	\
	modular-float.h		\
	modular-double.h	\
	modular-integer.h	\
	modular-inttype.h	\
	modular-log16.h		\
	modular-ruint.h		\
	modular-balanced.h 			\
	modular-balanced-int32.h 	\
	modular-balanced-int64.h 	\
	modular-balanced-float.h 	\
	modular-balanced-double.h 	\
	unparametric-operations.h	\
	unparametric.h		\
	zring.h				\
	montgomery.h	    \
	montgomery-int32.h	\
	montgomery-ruint.h	\
	modular-general.inl	\
	modular-int8.inl	\
	modular-int16.inl	\
	modular-int32.inl	\
	modular-int64.inl	\
	modular-uint8.inl	\
	modular-uint16.inl	\
	modular-uint32.inl	\
	modular-uint64.inl	\
	modular-float.inl	\
	modular-double.inl	\
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
	ring-interface.h

EXTRA_DIST=ring.doxy
