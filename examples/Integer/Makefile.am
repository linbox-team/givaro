all: $(EXTRA_PROGRAMS)
INCLUDES=-I$(top_srcdir)

AM_CPPFLAGS = -O7 -funroll-all-loops -felide-constructors -fstrict-aliasing -frerun-loop-opt -fexpensive-optimizations -fforce-addr -fforce-mem -fstrength-reduce -ffast-math $(GMP_CFLAGS) 
AM_CPPFLAGS += -I$(top_srcdir)/src/kernel/system -I$(top_srcdir)/src/kernel/memory -I$(top_srcdir)/src/kernel/integer -I$(top_srcdir)/src/kernel -I$(top_srcdir)/src/kernel/zpz -I$(top_srcdir)/src/kernel/bstruct

#LDADD = -L${libdir} -lgivaro $(GMP_LIBS)

LDADD = -L${top_srcdir}/src -lgivaro $(GMP_LIBS)

noinst_PROGRAMS=ifactor ifactor_lenstra igcd igcdext ilcm isproot nb_primes isprime nextprime order prevprime RSA_breaking RSA_keys_generator primitiveroot primitiveelement phi zeta zeta_inv iexponentiation RSA_decipher  RSA_encipher



ifactor_SOURCES = ifactor.C
ifactor_lenstra_SOURCES = ifactor_lenstra.C
igcd_SOURCES = igcd.C
zeta_SOURCES = zeta.C
zeta_inv_SOURCES = zeta_inv.C
phi_SOURCES = phi.C
igcdext_SOURCES = igcdext.C
ilcm_SOURCES = ilcm.C
isproot_SOURCES = isproot.C
primitiveroot_SOURCES = primitiveroot.C
primitiveelement_SOURCES = primitiveelement.C
nb_primes_SOURCES = nb_primes.C
isprime_SOURCES = isprime.C
nextprime_SOURCES = nextprime.C
order_SOURCES = order.C
prevprime_SOURCES = prevprime.C
RSA_breaking_SOURCES = RSA_breaking.C
RSA_keys_generator_SOURCES = RSA_keys_generator.C
RSA_encipher_SOURCES = RSA_encipher.C
RSA_decipher_SOURCES = RSA_decipher.C
iexponentiation_SOURCES = iexponentiation.C
