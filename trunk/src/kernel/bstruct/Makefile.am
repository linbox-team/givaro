INCLUDES=-I$(top_srcdir)

AM_CPPFLAGS= -I$(top_srcdir)/src/kernel/system -I$(top_srcdir)/src/kernel/memory

pkginclude_HEADERS=\
	givarray0.h		\
	givarrayallocator.h	\
	givarrayfixed.h		\
	givbits.h		\
	givelem.h		\
	givhashtable.h		\
	givlist0.h		\
	givstack.h		\
	givarray0.inl		\
	givbits.inl		\
	givhashtable.inl	\
	givlist0.inl		\
	givstack.inl		

noinst_LTLIBRARIES=libgivbstruct.la

libgivbstruct_la_SOURCES=givbits.C