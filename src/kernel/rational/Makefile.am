INCLUDES=-I$(top_srcdir)  ${GMP_VERSION}

AM_CPPFLAGS= -I$(top_srcdir)/src/kernel/integer -I$(top_srcdir)/src/kernel -I$(top_srcdir)/src/kernel/system -I$(top_srcdir)/src/kernel/memory

pkginclude_HEADERS=\
	givrational.h	\
	givrational.inl

noinst_LTLIBRARIES=libgivrational.la

libgivrational_la_SOURCES=\
	givratcstor.C	\
	givratio.C	\
	givratmisc.C	\
	givratmuldiv.C	\
	givrataddsub.C	\
	givratcpy.C	\
	givratcompare.C \
	givratreconstruct.C
