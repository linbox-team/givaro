# Copyright(c)'1994-2009 by The Givaro group
# This file is part of Givaro.
# Givaro is governed by the CeCILL-B license under French law
# and abiding by the rules of distribution of free software.
# see the COPYRIGHT file for more details.

AM_CPPFLAGS=-I$(top_srcdir)  ${GMP_VERSION}

AM_CXXFLAGS = @DEFAULT_CFLAGS@
AM_CPPFLAGS+= -I$(top_srcdir)/src/kernel/memory  -I$(top_srcdir)/src/kernel/system -I$(top_srcdir)/src/kernel -I$(top_srcdir)/src/kernel/ring

pkginclude_HEADERS=         \
	givinteger.h		    \
	random-integer.h		    \
	givintfactor.h		    \
	givintprime.h		    \
	givintnumtheo.h		    \
	givintrsa.h	  	        \
	givintrns.h		        \
	givintprime.inl		    \
	givintnumtheo.inl	    \
	givintrsa.inl		    \
	givintfactor.inl 	    \
	givintsqrootmod.h       \
	givintsqrootmod.inl     \
	givintrns_convert.inl	\
	givintrns_cstor.inl     \
	givspyinteger.h

noinst_LTLIBRARIES=libgivinteger.la

EXTRA_DIST=integer.doxy

libgivinteger_la_SOURCES=givintprime.C

#EXTRA_DIST= integer.doxy
