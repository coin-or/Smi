# Copyright (C) 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.

include $(top_srcdir)/BuildTools/Makemain.inc

########################################################################
#                               libSmi                                 #
########################################################################

# Name of the library compiled in this directory.
lib_LTLIBRARIES = libSmi.la

# List all source files for this library, including headers
libSmi_la_SOURCES = \
	SmiConfig.h \
	SmiCoreCombineRule.cpp SmiCoreCombineRule.hpp \
	SmiDiscreteDistribution.cpp SmiDiscreteDistribution.hpp \
	SmiLinearData.hpp \
	SmiScenarioTree.hpp \
	SmiScnData.cpp SmiScnData.hpp \
	SmiScnModel.cpp SmiScnModel.hpp \
	SmiMessage.cpp SmiMessage.hpp \
	SmiSmpsIO.cpp SmiSmpsIO.hpp

# List all additionally required libraries
libSmi_la_LIBADD = $(SMILIB_LFLAGS)

# Flags for the linker and preprocessor
AM_LDFLAGS = $(LT_LDFLAGS)
AM_CPPFLAGS = $(SMILIB_CFLAGS)

########################################################################
#                Headers that need to be installed                     #
########################################################################

# Here list all the header files that are required by a user of the library,
# and that therefore should be installed in 'include/coin-or'
includecoindir = $(pkgincludedir)
includecoin_HEADERS = \
	SmiCoreCombineRule.hpp \
	SmiDiscreteDistribution.hpp \
	SmiLinearData.hpp \
	SmiScenarioTree.hpp \
	SmiScnData.hpp \
	SmiMessage.hpp \
	SmiScnModel.hpp \
	SmiQuadratic.hpp

install-exec-local:
	$(install_sh_DATA) config_smi.h $(DESTDIR)$(includecoindir)/SmiConfig.h

uninstall-local:
	rm -f $(DESTDIR)$(includecoindir)/SmiConfig.h
