# Copyright(c)'1994-2009 by The Givaro group
# This file is part of Givaro.
# Givaro is governed by the CeCILL-B license under French law
# and abiding by the rules of distribution of free software.
# see the COPYRIGHT file for more details.
all: $(EXTRA_PROGRAMS)
INCLUDES=-I$(top_srcdir)

#OPTFLAGS = -O7 -funroll-all-loops -felide-constructors -fstrict-aliasing
#OPTFLAGS+= -frerun-loop-opt -fexpensive-optimizations
#OPTFLAGS+= -fomit-frame-pointer
#OPTFLAGS+= -fprefetch-loop-arrays -floop-optimize
#OPTFLAGS+= -malign-double
#OPTFLAGS+= -falign-loops -falign-jumps -falign-functions -falign-labels
#OPTFLAGS+= -fschedule-insns2
#OPTFLAGS+= -fforce-addr -fforce-mem -fstrength-reduce
#OPTFLAGS+= -ffast-math

# for icc ...
#OPTFLAGS = -fast -Ob2  -ipo_obj -unroll
#OPTFLAGS+= -parallel -par_report1
# icc for itanium2
#OPTFLAGS+= -tpp2 -mcpu=itanium2


AM_CXXFLAGS = @DEFAULT_CFLAGS@
AM_CPPFLAGS = $(OPTFLAGS) $(GMP_CFLAGS) -I$(top_srcdir)/src/kernel/system -I$(top_srcdir)/src/kernel/memory -I$(top_srcdir)/src/kernel/integer -I$(top_srcdir)/src/kernel -I$(top_srcdir)/src/library/poly1 -I$(top_srcdir)/src/kernel/zpz -I$(top_srcdir)/src/kernel/bstruct -I$(top_srcdir)/src/kernel/rational -I$(top_srcdir)/src/library/tools

LDADD = -L${top_srcdir}/src -L${top_srcdir}/src/.libs -lgivaro $(GMP_LIBS)
AM_LDFLAGS=-static


EXTRA_PROGRAMS=isirred isprimitive trunc_arith pol_arith pol_eval pol_factor interpolate PolynomialCRT
CLEANFILES=$(EXTRA_PROGRAMS)

interpolate_SOURCES   = interpolate.C
isirred_SOURCES       = isirred.C
isprimitive_SOURCES   = isprimitive.C
pol_arith_SOURCES     = pol_arith.C
trunc_arith_SOURCES   = trunc_arith.C
pol_eval_SOURCES      = pol_eval.C
pol_factor_SOURCES    = pol_factor.C
PolynomialCRT_SOURCES = PolynomialCRT.C

# for compilation of new examples
GIVARO_BIN=@bindir@
#  %:%.C
#  $(CXX) $(OPTFLAGS) ${INCLUDES} $(AM_CPPFLAGS) $*.C -o $*  $(shell $(GIVARO_BIN)/givaro-config --cflags ) $(LDADD)  $(shell $(GIVARO_BIN)/givaro-config --libs)
.C:
	$(CXX) $(OPTFLAGS) ${INCLUDES} $(AM_CPPFLAGS) $*.C -o $*  `$(GIVARO_BIN)/givaro-config --cflags` $(LDADD)  `$(GIVARO_BIN)/givaro-config --libs`
