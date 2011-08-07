# Copyright (C) 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Common Public License.

## $Id$

# Author:  Andreas Waechter           IBM    2006-04-13

AUTOMAKE_OPTIONS = foreign

########################################################################
#                               libSmi                                 #
########################################################################

# Name of the library compiled in this directory.  We want it to be installed
# in the 'lib' directory
lib_LTLIBRARIES = libSmi.la

# List all source files for this library, including headers
libSmi_la_SOURCES = \
	SmiCoreCombineRule.cpp SmiCoreCombineRule.hpp \
	SmiDiscreteDistribution.cpp SmiDiscreteDistribution.hpp \
	SmiLinearData.hpp \
	SmiScenarioTree.hpp \
	SmiScnData.cpp SmiScnData.hpp \
	SmiScnModel.cpp SmiScnModel.hpp \
	SmiMessage.cpp SmiMessage.hpp \
	SmiSmpsIO.cpp SmiSmpsIO.hpp

# This is for libtool (on Windows)
libSmi_la_LDFLAGS = $(LT_LDFLAGS)

########################################################################
#                            Additional flags                          #
########################################################################

# Here list all include flags, relative to this "srcdir" directory.  This
# "cygpath" stuff is necessary to compile with native compilers on Windows.
# "top_srcdir" refers to the basic directory for the main package that is
# being compiled.
AM_CPPFLAGS = $(SMI_CFLAGS)

# This line is necessary to allow VPATH compilation
DEFAULT_INCLUDES = -I. -I`$(CYGPATH_W) $(srcdir)`

########################################################################
#                Headers that need to be installed                     #
########################################################################

# Here list all the header files that are required by a user of the library,
# and that therefore should be installed in 'include/coin'
includecoindir = $(includedir)/coin
includecoin_HEADERS = \
	SmiCoreCombineRule.hpp \
	SmiDiscreteDistribution.hpp \
	SmiLinearData.hpp \
	SmiScenarioTree.hpp \
	SmiScnData.hpp \
	SmiMessage.hpp \
	SmiScnModel.hpp

install-exec-local:
	$(install_sh_DATA) config_smi.h $(DESTDIR)$(includecoindir)/SmiConfig.h

uninstall-local:
	rm -f $(DESTDIR)$(includecoindir)/SmiConfig.h
	