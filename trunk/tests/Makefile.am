# Copyright(c)'1994-2009 by The Givaro group
# This file is part of Givaro.
# Givaro is governed by the CeCILL-B license under French law
# and abiding by the rules of distribution of free software. 
# see the COPYRIGHT file for more details.

check:$(BASIC_TESTS)

SUBDIRS =
check: $(BASE_TESTS)

INCLUDES=-I$(top_srcdir) 
#AM_LDFLAGS = 
AM_CXXFLAGS = -g -Wall
AM_CPPFLAGS = $(OPTFLAGS) $(GMP_CFLAGS) -I$(top_srcdir)/src/kernel/system -I$(top_srcdir)/src/kernel/memory -I$(top_srcdir)/src/kernel/zpz -I$(top_srcdir)/src/kernel/integer -I$(top_srcdir)/src/kernel -I$(top_srcdir)/src/library/poly1 -I$(top_srcdir)/src/kernel/bstruct -I$(top_srcdir)/src/library/tools

LDADD = $(OPTLINKS) -L${libdir} -L${top_srcdir}/src -L${top_srcdir}/src/.libs -lgivaro $(GMP_LIBS)

BASIC_TESTS =   \
		test-random  \
		test-ifactor

CLEAN_FILES =   \
			 $(BASIC_TESTS)

EXTRA_PROGRAMS =   \
				 $(BASIC_TESTS)

TESTS =   \
	   $(BASIC_TESTS)

test_random_SOURCES= test-random.C
test_ifacor_SOURCES= test-ifactor.C

# for compilation of new tests (?)
GIVARO_BIN=@bindir@
%:%.C
	$(CXX) $(OPTFLAGS) ${INCLUDES} $(AM_CPPFLAGS) $*.C -o $*  $(shell $(GIVARO_BIN)/givaro-config --cflags ) $(LDADD)  $(shell $(GIVARO_BIN)/givaro-config --libs)

