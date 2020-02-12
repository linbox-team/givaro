# Copyright(c)'1994-2009 by The Givaro group
# This file is part of Givaro.
# Givaro is governed by the CeCILL-B license under French law
# and abiding by the rules of distribution of free software.
# see the COPYRIGHT file for more details.

#  all all-local: $(none)

SUBDIRS=Rational Integer FiniteField Polynomial RecInt Matrix

.PHONY:$(SUBDIRS)

examples: $(SUBDIRS)

$(SUBDIRS):
	$(MAKE) -C $@ examples

CLEANFILES:
	$(SUBDIRS)

