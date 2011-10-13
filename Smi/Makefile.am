# Copyright (C) 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Common Public License.

## $Id$

# Author:  Andreas Waechter           IBM    2006-04-13

AUTOMAKE_OPTIONS = foreign

########################################################################
#                          Subdirectories                              #
########################################################################

SUBDIRS = src

# We don't want to compile the test subdirectory, unless the test target is
# specified.  But we need to list it as subdirectory to make sure that it is
# included in the tarball

if ALWAYS_FALSE
  SUBDIRS += test
endif

########################################################################
#             Additional files to be included in tarball               #
########################################################################

# Here we need include all files that are not mentioned in other Makefiles

EXTRA_DIST = \
        examples/Makefile.in \
        examples/README \
        examples/stoch.cpp \
        test/SmiTestData/app0110.stoch \
        test/SmiTestData/app0110.cor \
        test/SmiTestData/app0110.time \
        test/SmiTestData/app0110R.stoch \
        test/SmiTestData/app0110R.cor \
        test/SmiTestData/app0110R.time \
        test/SmiTestData/bug.stoch \
        test/SmiTestData/bug.cor \
        test/SmiTestData/bug.time 

########################################################################
#                           Extra Targets                              #
########################################################################

test: all
	cd test; $(MAKE) test

unitTest: test

clean-local:
	if test -r test/Makefile; then cd test; $(MAKE) clean; fi

distclean-local:
	if test -r test/Makefile; then cd test; $(MAKE) distclean; fi

.PHONY: test unitTest 

########################################################################
#                  Installation of the addlibs file                    #
########################################################################

pkgconfiglibdir = $(libdir)/pkgconfig
pkgconfiglib_DATA = smi.pc

addlibsdir = $(DESTDIR)$(datadir)/coin/doc/Smi

install-data-hook:
	@$(mkdir_p) "$(addlibsdir)"
if COIN_HAS_PKGCONFIG
	PKG_CONFIG_PATH=@COIN_PKG_CONFIG_PATH@ \
	$(PKG_CONFIG) --libs smi > $(addlibsdir)/smi_addlibs.txt
else
if COIN_CXX_IS_CL
	echo "-libpath:`$(CYGPATH_W) @abs_lib_dir@` libSmi.lib @SMI_LIBS_INSTALLED@" > $(addlibsdir)/smi_addlibs.txt
else
	echo -L@abs_lib_dir@ -lSmi @SMI_LIBS_INSTALLED@ > $(addlibsdir)/smi_addlibs.txt
endif
endif

uninstall-hook:
	rm -f $(addlibsdir)/smi_addlibs.txt

########################################################################
#                         Maintainer Stuff                             #
########################################################################

# Files that are generated and should be cleaned with make distclean
DISTCLEANFILES = 

include BuildTools/Makemain.inc
