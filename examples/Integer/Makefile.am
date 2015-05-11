# Copyright(c)'1994-2009 by The Givaro group
# This file is part of Givaro.
# Givaro is governed by the CeCILL-B license under French law
# and abiding by the rules of distribution of free software.
# see the COPYRIGHT file for more details.
examples: $(EXTRA_PROGRAMS)
AM_CPPFLAGS=-I$(top_srcdir)

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
AM_CPPFLAGS += $(OPTFLAGS) $(GMP_CFLAGS)
AM_CPPFLAGS += -I$(top_srcdir)/src/kernel/system -I$(top_srcdir)/src/kernel/memory -I$(top_srcdir)/src/kernel/integer -I$(top_srcdir)/src/kernel -I$(top_srcdir)/src/kernel/field -I$(top_srcdir)/src/kernel/ring -I$(top_srcdir)/src/kernel/bstruct -I$(top_srcdir)/src/library/poly1 -I$(top_srcdir)/src/library/tools

LDADD = -L${top_srcdir}/src -L${top_srcdir}/src/.libs -lgivaro $(GMP_LIBS) $(LDFLAGS)
AM_LDFLAGS=-static

EXTRA_PROGRAMS=ifactor ifactor_lenstra igcd igcdext ilcm isproot nb_primes isprime nextprime order prevprime RSA_breaking RSA_keys_generator primitiveroot primitiveelement phi lambda lambda_inv iexponentiation RSA_decipher  RSA_encipher  ispower probable_primroot ProbLucas ModularSquareRoot

CLEANFILES=$(EXTRA_PROGRAMS)

ModularSquareRoot_SOURCES  = ModularSquareRoot.C
ProbLucas_SOURCES          = ProbLucas.C
probable_primroot_SOURCES  = probable_primroot.C
ispower_SOURCES            = ispower.C
ifactor_SOURCES            = ifactor.C
ifactor_lenstra_SOURCES    = ifactor_lenstra.C
igcd_SOURCES               = igcd.C
lambda_SOURCES             = lambda.C
lambda_inv_SOURCES         = lambda_inv.C
phi_SOURCES                = phi.C
igcdext_SOURCES            = igcdext.C
ilcm_SOURCES               = ilcm.C
isproot_SOURCES            = isproot.C
primitiveroot_SOURCES      = primitiveroot.C
primitiveelement_SOURCES   = primitiveelement.C
nb_primes_SOURCES          = nb_primes.C
isprime_SOURCES            = isprime.C
nextprime_SOURCES          = nextprime.C
order_SOURCES              = order.C
prevprime_SOURCES          = prevprime.C
RSA_breaking_SOURCES       = RSA_breaking.C
RSA_keys_generator_SOURCES = RSA_keys_generator.C
RSA_encipher_SOURCES       = RSA_encipher.C
RSA_decipher_SOURCES       = RSA_decipher.C
iexponentiation_SOURCES    = iexponentiation.C

# for compilation of new examples
GIVARO_BIN=@bindir@
%:%.C
	$(CXX) $(OPTFLAGS) $(AM_CPPFLAGS) $(AM_CXXFLAGS) $(CXXFLAGS) $*.C -o $*  `$(GIVARO_BIN)/givaro-config --cflags` $(LDADD)  `$(GIVARO_BIN)/givaro-config --libs`
