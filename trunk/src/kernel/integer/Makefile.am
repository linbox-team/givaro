INCLUDES=-I$(top_srcdir)  ${GMP_VERSION}

AM_CPPFLAGS= -I$(top_srcdir)/src/kernel/memory  -I$(top_srcdir)/src/kernel/system -I$(top_srcdir)/src/kernel

pkginclude_HEADERS=\
	givinteger.h		\
	givintfactor.h		\
	givintprime.h		\
	givintnumtheo.h		\
	givintrsa.h		\
	givintrns.h		\
	givinteger.inl		\
	givintprime.inl		\
	givintnumtheo.inl	\
	givintrsa.inl		\
	givintfactor.inl 	\
        givintrns_convert.inl	\
	givintrns_cstor.inl

noinst_LTLIBRARIES=libgivinteger.la

libgivinteger_la_SOURCES=givintprime.C
