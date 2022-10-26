# Copyright(c)'1994-2009 by The Givaro group
# This file is part of Givaro.
# Givaro is governed by the CeCILL-B license under French law
# and abiding by the rules of distribution of free software.
# see the COPYRIGHT file for more details.

AM_CPPFLAGS=-I$(top_srcdir)

AM_CXXFLAGS = $(GIVARO_CXXFLAGS)

AM_LDFLAGS=$(LDFLAGS)

SUBDIRS=kernel library

lib_LTLIBRARIES= libgivaro.la

libgivaro_la_SOURCES = dummy.C

libgivaro_la_LIBADD= \
	kernel/gmp++/libgmppp.la kernel/bstruct/libgivbstruct.la kernel/integer/libgivinteger.la kernel/memory/libgivmemory.la kernel/rational/libgivrational.la kernel/system/libgivsystem.la kernel/field/libgivfield.la library/tools/libgivtools.la library/poly1/libgivpoly1.la $(GMP_LIBS) 

# soname of libgivaro
libgivaro_la_LDFLAGS = -version-info 11:1:2
libgivaro_la_LDFLAGS+= -no-undefined
libgivaro_la_LDFLAGS+=$(LDFLAGS)
