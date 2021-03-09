# Copyright(c)'1994-2009 by The Givaro group
# This file is part of Givaro.
# Givaro is governed by the CeCILL-B license under French law
# and abiding by the rules of distribution of free software.
# see the COPYRIGHT file for more details.
examples: $(EXTRA_PROGRAMS)
AM_CPPFLAGS=-I$(top_builddir)

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


AM_CXXFLAGS = $(GIVARO_CXXFLAGS)
AM_CPPFLAGS += $(OPTFLAGS) -I$(top_builddir)/src/kernel/system -I$(top_builddir)/src/kernel/memory -I$(top_builddir)/src/kernel/integer -I$(top_builddir)/src/kernel -I$(top_builddir)/src/library/poly1 -I$(top_builddir)/src/kernel/field -I$(top_builddir)/src/kernel/ring -I$(top_builddir)/src/kernel/bstruct -I$(top_builddir)/src/kernel/rational -I$(top_builddir)/src/library/tools -I$(top_srcdir)/src/kernel

LDADD = $(top_builddir)/src/libgivaro.la $(GMP_LIBS) $(LDFLAGS)
AM_LDFLAGS=-static


EXTRA_PROGRAMS=isirred isprimitive trunc_arith pol_arith pol_eval pol_factor interpolate PolynomialCRT highorder bivariate isgenerator AES_SubBytes

CLEANFILES=$(EXTRA_PROGRAMS)

interpolate_SOURCES   = interpolate.C
isirred_SOURCES       = isirred.C
isgenerator_SOURCES       = isgenerator.C
isprimitive_SOURCES   = isprimitive.C
pol_arith_SOURCES     = pol_arith.C
trunc_arith_SOURCES   = trunc_arith.C
pol_eval_SOURCES      = pol_eval.C
pol_factor_SOURCES    = pol_factor.C
PolynomialCRT_SOURCES = PolynomialCRT.C
highorder_SOURCES     = highorder.C
bivariate_SOURCES     = bivariate.C
AES_SubBytes_SOURCES     = AES_SubBytes.C

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
