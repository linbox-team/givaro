# Copyright(c)'1994-2009 by The Givaro group
# This file is part of Givaro.
# Givaro is governed by the CeCILL-B license under French law
# and abiding by the rules of distribution of free software.
# see the COPYRIGHT file for more details.

ACLOCAL_AMFLAGS = -I macros

SUBDIRS=src macros tests docs

bin_SCRIPTS=givaro-config givaro-makefile

include_HEADERS=givaro-config.h

EXTRA_DIST=COPYRIGHT Licence_CeCILL-B_V1-en.txt Licence_CeCILL-B_V1-fr.txt \
		   givaro.doxy

DISTCLEANFILES = src/kernel/system/givaro
DISTCLEANFILES += src/kernel/bstruct/givaro
DISTCLEANFILES += src/kernel/integer/givaro
DISTCLEANFILES += src/kernel/memory/givaro
DISTCLEANFILES += src/kernel/rational/givaro
DISTCLEANFILES += src/kernel/zpz/givaro
DISTCLEANFILES += src/kernel/gmp++/gmp++
DISTCLEANFILES += src/library/matrix/givaro
DISTCLEANFILES += src/library/poly1/givaro
DISTCLEANFILES += src/library/tools/givaro
DISTCLEANFILES += src/library/vector/givaro

.PHONY: examples

examples:
	$(MAKE) && cd examples && $(MAKE) $(AM_MAKEFLAGS) examples

docs:docs/givaro-html/index.html

docs/givaro-html/index.html:
	(cd docs; make docs)

#docs_dev:docs/givaro-dev-html/index.html

#doc/givaro-dev-html/index.html:
#(cd docs; make docs_dev)



