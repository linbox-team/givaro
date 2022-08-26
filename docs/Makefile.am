# Copyright(c)'2010 by The Givaro group
# This file is part of Givaro.
# Givaro is governed by the CeCILL-B license under French law
# and abiding by the rules of distribution of free software.
# see the COPYRIGHT file for more details.


if GIVARO_BUILD_DOC
USE_TARGETS = docs
INSTALL_TARGETS = install-doc
endif

docdir=$(GIVARO_DOC_PATH)

#man1_MANS = givaro-config.1

all all-local: $(USE_TARGETS)

install install-data-local: $(USE_TARGETS) $(INSTALL_TARGETS)

if GIVARO_DOXYGEN_FOUND

docs:
	sed -i 's/^\\version.*/\\version\ $(VERSION)/' mainpage.doxy
	if test -d givaro-html ; then echo exists; else mkdir givaro-html ; fi
	if test -d givaro-dev-html ; then echo exists; else mkdir givaro-dev-html ; fi
	cp ../INSTALL givaro-html/
	cp ../COPYING givaro-html/
	cp ../AUTHORS givaro-html/
	doxygen Doxyfile

#  if test -d givaro-dev-html ; then echo exists; else mkdir givaro-dev-html ; fi
#  cp index-dev.html givaro-dev-html/index.html
#  cp tutorial.html givaro-html/
#  cp install-dist.html givaro-html/
#  cp install-dev.html givaro-html/

docs_dev:
	make docs
	cp ../INSTALL givaro-dev-html/
	cp ../COPYING givaro-dev-html/
	cp ../AUTHORS givaro-dev-html/
	doxygen DoxyfileDev

#  cp tutorial.html givaro-dev-html/
#  cp install-dist.html givaro-dev-html/
#  cp install-dev.html givaro-dev-html/

install-doc:
	mkdir -p               $(DESTDIR)/$(docdir)
	cp -rp givaro-html     $(DESTDIR)/$(docdir)/givaro-html
	cp -rp givaro-dev-html $(DESTDIR)/$(docdir)/givaro-dev-html
	cp -p givaro.html      $(DESTDIR)/$(docdir)/givaro.html

uninstall-hook:
	(test -d "$(DESTDIR)/$(docdir)/givaro-html" && rm -r "$(DESTDIR)/$(docdir)/givaro-html") || true
	(test -d "$(DESTDIR)/$(docdir)/givaro-dev-html" && rm -r "$(DESTDIR)/$(docdir)/givaro-dev-html") || true
	(test -f "$(DESTDIR)/$(docdir)/givaro.html" && rm "$(DESTDIR)/$(docdir)/givaro.html") || true

else

docs:
	@echo
	@echo "*** ERROR ***"
	@echo "you need doxygen to build documentation"
	@echo

endif

EXTRA_DIST= \
			Doxyfile.in    \
			DoxyfileDev.in \
			mainpage.doxy  \
			givaro.html
# \
#			doc.doxy            \
#			tutorial.doxy       \
#			givaro.html         \
#			givaro-config.1 \
#tutorial.html       \
#	install-dev.html    \
#	index-dev.html      \
#	install-dist.html


clean-local :
	rm -rf givaro-html
	rm -rf givaro-dev-html


