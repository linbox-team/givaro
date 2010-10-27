# Copyright(c)'1994-2009 by The Givaro group
# This file is part of Givaro.
# Givaro is governed by the CeCILL-B license under French law
# and abiding by the rules of distribution of free software. 
# see the COPYRIGHT file for more details.
all: $(EXTRA_PROGRAMS)

INCLUDES=-I$(top_srcdir)

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

AM_CPPFLAGS = $(OPTFLAGS) $(GMP_CFLAGS) -I$(top_srcdir)/src/kernel/system -I$(top_srcdir)/src/kernel/memory -I$(top_srcdir)/src/kernel/zpz -I$(top_srcdir)/src/kernel/integer -I$(top_srcdir)/src/kernel -I$(top_srcdir)/src/library/poly1 -I$(top_srcdir)/src/kernel/bstruct -I$(top_srcdir)/src/library/tools

LDADD = $(OPTLINKS) -L${libdir} -L${top_srcdir}/src -L${top_srcdir}/src/.libs -lgivaro $(GMP_LIBS)


EXTRA_PROGRAMS=gfq_atomic Test_Extension zpz_atomic ff_arith exponentiation domain_to_operatorstyle all_field ChineseRemaindering GF128 GFirreducible

GFirreducible_SOURCES=GFirreducible.C
GF128_SOURCES=GF128.C
ChineseRemaindering_SOURCES=ChineseRemaindering.C
domain_to_operatorstyle_SOURCES=domain_to_operatorstyle.C
all_field_SOURCES=all_field.C
gfq_atomic_SOURCES=gfq_atomic.C
Test_Extension_SOURCES=Test_Extension.C
zpz_atomic_SOURCES=zpz_atomic.C
ff_arith_SOURCES=ff_arith.C
exponentiation_SOURCES=exponentiation.C

# for compilation of new examples
GIVARO_BIN=@bindir@
%:%.C
	$(CXX) $(OPTFLAGS) ${INCLUDES} $(AM_CPPFLAGS) $*.C -o $*  $(shell $(GIVARO_BIN)/givaro-config --cflags ) $(LDADD)  $(shell $(GIVARO_BIN)/givaro-config --libs)
