# Copyright(c)'1994-2009 by The Givaro group
# This file is part of Givaro.
# Givaro is governed by the CeCILL-B license under French law
# and abiding by the rules of distribution of free software.
# see the COPYRIGHT file for more details.

ACLOCAL_AMFLAGS = -I macros

SUBDIRS=src macros tests docs

include_HEADERS=givaro-config.h

docs:docs/givaro-html/index.html

docs/givaro-html/index.html:
	(cd docs; make docs)

.PHONY:examples

examples:
	(cd examples; make examples)


bin_SCRIPTS=givaro-config givaro-makefile

EXTRA_DIST=COPYRIGHT Licence_CeCILL-B_V1-en.txt Licence_CeCILL-B_V1-fr.txt \
		   givaro.doxy

DIST_SUBDIRS=examples ${SUBDIRS}

DISTCLEANFILES = src/kernel/system/givaro   \
				 src/kernel/bstruct/givaro  \
				 src/kernel/integer/givaro  \
				 src/kernel/memory/givaro   \
				 src/kernel/rational/givaro \
				 src/kernel/zpz/givaro      \
				 src/kernel/gmp++/gmp++     \
				 src/library/matrix/givaro  \
				 src/library/poly1/givaro   \
				 src/library/tools/givaro   \
				 src/library/vector/givaro


# docs_dev:docs/givaro-dev-html/index.html

#  doc/givaro-dev-html/index.html:
#  (cd docs; make docs_dev)

VERSION=3.5.0

