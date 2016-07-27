# Copyright(c)'1994-2009 by The Givaro group
# This file is part of Givaro.
# Givaro is governed by the CeCILL-B license under French law
# and abiding by the rules of distribution of free software.
# see the COPYRIGHT file for more details.
examples: $(EXTRA_PROGRAMS)
AM_CPPFLAGS=-I$(top_srcdir)
AM_CXXFLAGS = @DEFAULT_CFLAGS@
AM_CPPFLAGS += $(OPTFLAGS) $(GMP_CFLAGS)
AM_CPPFLAGS += -I$(top_srcdir)/src/kernel

#OPTFLAGS=
#OPTLINKS=

LDADD = $(OPTLINKS) -L${top_srcdir}/src -L${top_srcdir}/src/.libs -lgivaro $(GMP_LIBS) $(LDFLAGS)
AM_LDFLAGS=-static


EXTRA_PROGRAMS=recint-iterator rsa extended-int-types
CLEANFILES=$(EXTRA_PROGRAMS)

recint_iterator_SOURCES  = recint-iterator.C
extended_int_types_SOURCES  = extended-int-types.C
rsa_SOURCES       = rsa.C

# for compilation of new examples
GIVARO_BIN=@bindir@
examp_comp = $(CXX) $(OPTFLAGS) $(AM_CPPFLAGS) $(AM_CXXFLAGS) $< -o $*  $(shell $(GIVARO_BIN)/givaro-config --cflags ) $(LDADD) $(shell $(GIVARO_BIN)/givaro-config --libs)

%:%.C
	$(examp_comp)

%:%.cpp
	$(examp_comp)

