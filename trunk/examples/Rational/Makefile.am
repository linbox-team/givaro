check: $(EXTRA_PROGRAMS)
INCLUDES=-I$(top_srcdir)

AM_CPPFLAGS = -O7 -funroll-all-loops -felide-constructors -fstrict-aliasing -frerun-loop-opt -fexpensive-optimizations -fforce-addr -fforce-mem -fstrength-reduce -ffast-math $(GMP_CFLAGS) 
AM_CPPFLAGS += -I$(top_srcdir)/src/kernel/system -I$(top_srcdir)/src/kernel/memory -I$(top_srcdir)/src/kernel/integer -I$(top_srcdir)/src/kernel -I$(top_srcdir)/src/kernel/zpz -I$(top_srcdir)/src/kernel/bstruct -I$(top_srcdir)/src/kernel/rational

#LDADD = -L${libdir} -lgivaro $(GMP_LIBS)

LDADD = -L${top_srcdir}/src -lgivaro $(GMP_LIBS)

EXTRA_PROGRAMS=iratrecon

iratrecon_SOURCES = iratrecon.C
