INCLUDES=-I$(top_srcdir)

AM_CPPFLAGS= -I$(top_srcdir)/src/kernel/bstruct -I$(top_srcdir)/src/kernel/system -I$(top_srcdir)/src/kernel/memory

pkginclude_HEADERS=\
	givindeter.h		\
	givpoly1.h		\
	givdegree.h		\
	givinterp.h		\
	givpoly1dense.h		\
	givpoly1factor.h	\
	givpoly1padic.h		\
	givpoly1denseops.inl	\
	givpoly1addsub.inl	\
	givpoly1muldiv.inl	\
	givpoly1axpy.inl	\
	givpoly1gcd.inl		\
	givpoly1kara.inl	\
	givpoly1cstor.inl	\
	givpoly1io.inl		\
	givpoly1misc.inl	\
	givpoly1sqrfree.inl	\
	givpoly1cyclo.inl	\
	givpoly1factor.inl	\
	givpoly1proot.inl

noinst_LTLIBRARIES=libgivpoly1.la

libgivpoly1_la_SOURCES=givindeter.C givdegree.C
