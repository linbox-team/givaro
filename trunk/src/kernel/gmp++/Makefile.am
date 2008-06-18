INCLUDES=-I$(top_srcdir) ${GMP_VERSION}

pkgincludesubdir=$(includedir)/gmp++

AM_CPPFLAGS= -I$(top_srcdir)/src/kernel/memory -I$(top_srcdir)/src/kernel/system

pkgincludesub_HEADERS=\
	gmp++.h		\
	gmp++_int.h	\
	gmp++_int.inl	\
	gmp++.C		\
	gmp++_int.C

lib_LTLIBRARIES=libgmpxx.la

libgmpxx_la_SOURCES=\
	gmp++_int_div.C		\
	gmp++_int_mod.C		\
	gmp++_int_sub.C		\
	gmp++_int_compare.C	\
	gmp++_int_gcd.C		\
	gmp++_int_io.C		\
	gmp++_int_mul.C		\
	gmp++_int_add.C		\
	gmp++_int_cstor.C	\
	gmp++_int_misc.C	\
	gmp++_int_pow.C
