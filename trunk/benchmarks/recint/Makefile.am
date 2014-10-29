# Copyright(c)'1994-2009 by The Givaro group
# This file is part of Givaro.
# Givaro is governed by the CeCILL-B license under French law
# and abiding by the rules of distribution of free software.
# see the COPYRIGHT file for more details.

benchmarks: $(EXTRA_PROGRAMS)

AM_CPPFLAGS=-I$(top_srcdir)

OPTFLAGS=
OPTLINKS=

AM_CXXFLAGS = @DEFAULT_CFLAGS@
AM_CPPFLAGS += $(OPTFLAGS) $(GMP_CFLAGS) -I$(top_srcdir)/src/kernel/system -I$(top_srcdir)/src/kernel/recint

LDADD = $(OPTLINKS) -L${top_srcdir}/src -L${top_srcdir}/src/.libs -lgivaro $(GMP_LIBS) $(LDFLAGS)
AM_LDFLAGS=-static

EXTRA_PROGRAMS=exp inv_arazi inv_naive mul
CLEANFILES=$(EXTRA_PROGRAMS)

exp_SOURCES = exp.C
inv_arazi_SOURCES = inv_arazi.C
inv_naive_SOURCES = inv_naive.C
mul_SOURCES = mul.C

