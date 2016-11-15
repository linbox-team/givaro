# Copyright(c)'1994-2009 by The Givaro group
# This file is part of Givaro.
# Givaro is governed by the CeCILL-B license under French law
# and abiding by the rules of distribution of free software.
# see the COPYRIGHT file for more details.
examples: $(EXTRA_PROGRAMS)

AM_CPPFLAGS=-I$(top_srcdir)

OPTFLAGS=
OPTLINKS=

# for gcc ...
#OPTFLAGS+= -O7 -funroll-all-loops -felide-constructors -fstrict-aliasing
#OPTFLAGS+= -frerun-loop-opt -fexpensive-optimizations
#OPTFLAGS+= -fomit-frame-pointer
#OPTFLAGS+= -fprefetch-loop-arrays -floop-optimize
#OPTFLAGS+= -malign-double
#OPTFLAGS+= -falign-loops -falign-jumps -falign-functions -falign-labels
#OPTFLAGS+= -fschedule-insns2
#OPTFLAGS+= -fforce-addr -fforce-mem -fstrength-reduce
#OPTFLAGS+= -ffast-math

# for icc ...
#OPTFLAGS+= -fast -Ob2  -ipo_obj -unroll
#OPTFLAGS+= -parallel -par_report1
#OPTLINKS+= -parallel
# icc for itanium2
#OPTFLAGS+= -tpp2 -mcpu=itanium2

AM_CXXFLAGS = @DEFAULT_CFLAGS@
AM_CPPFLAGS += $(OPTFLAGS) $(GMP_CFLAGS) -I$(top_srcdir)/src/kernel/system -I$(top_srcdir)/src/kernel/memory -I$(top_srcdir)/src/kernel/field -I$(top_srcdir)/src/kernel/ring -I$(top_srcdir)/src/kernel/integer -I$(top_srcdir)/src/kernel -I$(top_srcdir)/src/library/poly1 -I$(top_srcdir)/src/kernel/bstruct -I$(top_srcdir)/src/library/tools

LDADD = $(OPTLINKS) $(top_srcdir)/src/libgivaro.la $(GMP_LIBS) $(LDFLAGS)
AM_LDFLAGS=-static


EXTRA_PROGRAMS=gfq_atomic GaloisFields Test_Extension zpz_atomic ff_arith exponentiation domain_to_operatorstyle all_field GF128 GFirreducible

CLEANFILES=$(EXTRA_PROGRAMS)

GFirreducible_SOURCES           = GFirreducible.C
GF128_SOURCES                   = GF128.C
domain_to_operatorstyle_SOURCES = domain_to_operatorstyle.C
all_field_SOURCES               = all_field.C
gfq_atomic_SOURCES              = gfq_atomic.C
GaloisFields_SOURCES             = GaloisFields.C
Test_Extension_SOURCES          = Test_Extension.C
zpz_atomic_SOURCES              = zpz_atomic.C
ff_arith_SOURCES                = ff_arith.C
exponentiation_SOURCES          = exponentiation.C

# for compilation of new examples
GIVARO_BIN=@bindir@
examp_comp = $(CXX) $(OPTFLAGS) $(AM_CPPFLAGS) $(AM_CXXFLAGS) $(CXXFLAGS) $*.C -o $*  `$(GIVARO_BIN)/givaro-config --cflags` $(LDADD)  `$(GIVARO_BIN)/givaro-config --libs`

%:%.C
	$(examp_comp)

%:%.cpp
	$(examp_comp)

