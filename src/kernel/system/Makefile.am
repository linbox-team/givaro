INCLUDES=-I$(top_srcdir)

AM_CPPFLAGS= -I$(top_srcdir)/src/kernel/memory

pkginclude_HEADERS=\
	givbasictype.h	\
	givconfig.h	\
	giverror.h	\
	givgenarith.h	\
	givinit.h	\
	givmodule.h	\
	givperf.h	\
	givpower.h	\
	givrandom.h	\
	givtimer.h

noinst_LTLIBRARIES=libgivsystem.la

libgivsystem_la_SOURCES=\
	givbasictype.C	\
	giverror.C	\
	givinit.C	\
	givmodule.C	\
	givtimer.C
