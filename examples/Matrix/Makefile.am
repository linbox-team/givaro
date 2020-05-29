# Copyright(c)'1994-2009 by The Givaro group
# This file is part of Givaro.
# Givaro is governed by the CeCILL-B license under French law
# and abiding by the rules of distribution of free software.
# see the COPYRIGHT file for more details.
examples: $(EXTRA_PROGRAMS)
AM_CPPFLAGS =-I$(top_builddir)
AM_CXXFLAGS = $(GIVARO_CXXFLAGS)
AM_CPPFLAGS += $(OPTFLAGS) -I$(top_builddir)/src/kernel/system -I$(top_builddir)/src/kernel/memory -I$(top_builddir)/src/kernel/integer -I$(top_builddir)/src/kernel -I$(top_builddir)/src/library/tools -I$(top_builddir)/src/library/vector -I$(top_builddir)/src/library/poly1 -I$(top_builddir)/src/kernel/field -I$(top_builddir)/src/kernel/ring -I$(top_builddir)/src/kernel/bstruct -I$(top_builddir)/src/kernel/rational -I$(top_builddir)/src/library/matrix -I$(top_srcdir)/src/kernel

LDADD = $(top_builddir)/src/libgivaro.la $(GMP_LIBS) $(LDFLAGS)
AM_LDFLAGS=-static


EXTRA_PROGRAMS=f4n3

CLEANFILES=$(EXTRA_PROGRAMS)

f4n3_SOURCES   = f4n3.cpp

#external flags
AM_CPPFLAGS += $(OPTFLAGS)
LDADD += $(OPTLIBES)

# for compilation of new examples
define comp_new_examp
	$(AM_V_CXX)$(CXXCOMPILE) $(OPTFLAGS) -c -o $@.$(OBJEXT) $<
	$(AM_V_CXXLD)$(CXXLINK) $@.$(OBJEXT) $(LDADD)
endef

%:%.C
	$(comp_new_examp)

%:%.cpp
	$(comp_new_examp)
