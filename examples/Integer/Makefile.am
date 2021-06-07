# Copyright(c)'1994-2009 by The Givaro group
# This file is part of Givaro.
# Givaro is governed by the CeCILL-B license under French law
# and abiding by the rules of distribution of free software.
# see the COPYRIGHT file for more details.
examples: $(EXTRA_PROGRAMS)
AM_CPPFLAGS=-I$(top_builddir)

AM_CXXFLAGS = $(GIVARO_CXXFLAGS)
AM_CPPFLAGS += -I$(top_builddir)/src/kernel/system -I$(top_builddir)/src/kernel/memory -I$(top_builddir)/src/kernel/integer -I$(top_builddir)/src/kernel -I$(top_builddir)/src/kernel/field -I$(top_builddir)/src/kernel/ring -I$(top_builddir)/src/kernel/bstruct -I$(top_builddir)/src/library/poly1 -I$(top_builddir)/src/library/tools -I$(top_srcdir)/src/kernel

LDADD = $(top_builddir)/src/libgivaro.la $(GMP_LIBS) $(LDFLAGS)
AM_LDFLAGS=-static

EXTRA_PROGRAMS=ifactor ifactor_lenstra igcd igcdext ilcm isproot nb_primes isprime nextprime order prevprime RSA_breaking RSA_keys_generator primitiveroot primitiveelement phi lambda lambda_inv iexponentiation RSA_decipher  RSA_encipher  ispower probable_primroot ProbLucas ModularSquareRoot Fibonacci

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
Fibonacci_SOURCES    = Fibonacci.C

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
