check: $(EXTRA_PROGRAMS)
INCLUDES=-I$(top_srcdir)

AM_CPPFLAGS = -O7 -funroll-all-loops -felide-constructors -fstrict-aliasing -frerun-loop-opt -fexpensive-optimizations -fforce-addr -fforce-mem -fstrength-reduce -ffast-math $(GMP_CFLAGS) -I$(top_srcdir)/src/kernel/system -I$(top_srcdir)/src/kernel/memory -I$(top_srcdir)/src/kernel/integer -I$(top_srcdir)/src/kernel -I$(top_srcdir)/src/library/poly1 -I$(top_srcdir)/src/kernel/zpz -I$(top_srcdir)/src/kernel/bstruct -I$(top_srcdir)/src/library/tools

LDADD = -L${libdir} -L${top_srcdir}/src -L${top_srcdir}/src/.libs -lgivaro $(GMP_LIBS)


EXTRA_PROGRAMS=isirred pol_arith pol_eval

isirred_SOURCES = isirred.C
pol_arith_SOURCES = pol_arith.C
pol_eval_SOURCES = pol_eval.C
