# Copyright(c)'1994-2009 by The Givaro group
# This file is part of Givaro.
# Givaro is governed by the CeCILL-B license under French law
# and abiding by the rules of distribution of free software.
# see the COPYRIGHT file for more details.

SUBDIRS =
check: $(BASE_TESTS)

INCLUDES=-I$(top_srcdir)
#AM_LDFLAGS =
#  AM_CXXFLAGS = -g -Wall -O0
AM_CXXFLAGS = @TESTS_CFLAGS@
AM_CPPFLAGS = $(OPTFLAGS) $(GMP_CFLAGS) -I$(top_srcdir)/src/kernel/system -I$(top_srcdir)/src/kernel/memory -I$(top_srcdir)/src/kernel/zpz -I$(top_srcdir)/src/kernel/integer -I$(top_srcdir)/src/kernel -I$(top_srcdir)/src/library/poly1 -I$(top_srcdir)/src/kernel/bstruct -I$(top_srcdir)/src/library/tools

LDADD = $(OPTLINKS) -L${top_srcdir}/src -L${top_srcdir}/src/.libs -lgivaro $(GMP_LIBS)
AM_LDFLAGS=-static

BASIC_TESTS =          \
		test-random    \
		test-ifactor   \
		test-ffarith   \
		test-ringarith \
		test-mod       \
		test-modsqroot \
		test-trunc     \
		test-crt       \
		test-geom      \
		test-integer

CLEANFILES =           \
		$(BASIC_TESTS)

EXTRA_PROGRAMS =       \
		$(BASIC_TESTS)

TESTS =                \
	    $(BASIC_TESTS)


test_random_SOURCES   = test-random.C
test_ifactor_SOURCES  = test-ifactor.C
test_ffarith_SOURCES  = test-ffarith.C
test_ringarith_SOURCES= test-ringarith.C
test_mod_SOURCES      = test-mod.C
test_modsqroot_SOURCES= test-modsqroot.C
test_trunc_SOURCES    = test-trunc.C
test_crt_SOURCES      = test-crt.C
test_geom_SOURCES     = test-geom.C
test_integer_SOURCES  = test-integer.C


# for compilation of new tests
GIVARO_BIN=@bindir@
#  %:%.C
#  $(CXX) $(OPTFLAGS) ${INCLUDES} $(AM_CPPFLAGS) $*.C -o $*  $(shell $(GIVARO_BIN)/givaro-config --cflags ) $(LDADD)  $(shell $(GIVARO_BIN)/givaro-config --libs)
.C:
	$(CXX) $(OPTFLAGS) ${INCLUDES} $(AM_CPPFLAGS) $*.C -o $*  `$(GIVARO_BIN)/givaro-config --cflags` $(LDADD)  `$(GIVARO_BIN)/givaro-config --libs`
