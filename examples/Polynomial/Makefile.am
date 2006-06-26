check: $(EXTRA_PROGRAMS)
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


AM_CPPFLAGS = $(OPTFLAGS) $(GMP_CFLAGS) -I$(top_srcdir)/src/kernel/system -I$(top_srcdir)/src/kernel/memory -I$(top_srcdir)/src/kernel/integer -I$(top_srcdir)/src/kernel -I$(top_srcdir)/src/library/poly1 -I$(top_srcdir)/src/kernel/zpz -I$(top_srcdir)/src/kernel/bstruct -I$(top_srcdir)/src/library/tools

LDADD = -L${libdir} -L${top_srcdir}/src -L${top_srcdir}/src/.libs -lgivaro $(GMP_LIBS)


EXTRA_PROGRAMS=isirred isprimitive pol_arith pol_eval pol_factor interpolate

interpolate_SOURCES = interpolate.C
isirred_SOURCES = isirred.C
isprimitive_SOURCES = isprimitive.C
pol_arith_SOURCES = pol_arith.C
pol_eval_SOURCES = pol_eval.C
pol_factor_SOURCES = pol_factor.C

# for compilation of new examples
GIVARO_BIN=@bindir@
%:%.C 
	$(CXX) $(OPTFLAGS) ${INCLUDES} $(AM_CPPFLAGS) $*.C -o $*  $(shell $(GIVARO_BIN)/givaro-config --cflags ) $(LDADD)  $(shell $(GIVARO_BIN)/givaro-config --libs)
