# Copyright(c)'1994-2009 by The Givaro group
# This file is part of Givaro.
# Givaro is governed by the CeCILL-B license under French law
# and abiding by the rules of distribution of free software.
# see the COPYRIGHT file for more details.

SUBDIRS =
benchmarks: $(EXTRA_PROGRAMS)

AM_CPPFLAGS=-I$(top_builddir)

PERFPUBLISHERFILE=benchmarks-report.xml

OPTFLAGS=
OPTLINKS=

AM_CXXFLAGS = $(GIVARO_CXXFLAGS) $(OPTFLAGS)
AM_CPPFLAGS += -I$(top_builddir)/src/kernel/system -I$(top_builddir)/src/kernel/recint -I$(top_builddir)/src/kernel/integer -I$(top_builddir)/src/kernel/gmp++ -I$(top_builddir)/src/kernel/ring

LDADD = $(OPTLINKS) $(top_builddir)/src/libgivaro.la $(GMP_LIBS)
AM_LDFLAGS=-static

BENCHMARKS=				\
	benchmark-recint_exp			\
	benchmark-recint_inv_arazi	\
	benchmark-recint_inv_naive	\
	benchmark-recint_mul

CLEANFILES = $(EXTRA_PROGRAMS) $(PERFPUBLISHERFILE)

EXTRA_PROGRAMS = $(BENCHMARKS)

benchmark_recint_exp_SOURCES = benchmark-recint_exp.C
benchmark_recint_inv_arazi_SOURCES = benchmark-recint_inv_arazi.C
benchmark_recint_inv_naive_SOURCES = benchmark-recint_inv_naive.C
benchmark_recint_mul_SOURCES = benchmark-recint_mul.C

# Perfpublisher script interaction - AB 2014/11/17
perfpublisher:
	+./perfpublisher.sh "$(PERFPUBLISHERFILE)" "$(BENCHMARKS)" "$(CXX)"
	mv $(PERFPUBLISHERFILE) ..
	make clean
	mv ../$(PERFPUBLISHERFILE) .

#=====
# For compilation of new benchmarks

GIVARO_BIN=@bindir@
%:%.C
	$(CXX) $(OPTFLAGS) $(AM_CPPFLAGS) $(AM_CXXFLAGS) $(CXXFLAGS) $*.C -o $* `$(GIVARO_BIN)/givaro-config --cflags` $(LDADD) `$(GIVARO_BIN)/givaro-config --libs`

