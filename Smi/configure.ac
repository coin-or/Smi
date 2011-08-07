# Copyright (C) 2006 International Business Machines.
# All Rights Reserved.
# This file is distributed under the Common Public License.

## $Id$

# Author:  Andreas Waechter            IBM    2006-04-13

#############################################################################
#                       Names and other basic things                        #
#############################################################################

AC_PREREQ(2.59)

AC_INIT([Smi],[trunk],[smi@list.coin-or.org])

AC_COPYRIGHT([
Copyright 2006 International Business Machines and others.
All Rights Reserved.
This file is part of the open source package Smi which is distributed
under the Common Public License.])

# List one file in the package so that the configure script can test
# whether the package is actually there
AC_CONFIG_SRCDIR(src/SmiScnModel.hpp)

# Where should everything be installed by default?  Here, we want it
# to be installed directly in 'bin', 'lib', 'include' subdirectories
# of the directory where configure is run.  The default would be
# /usr/local.
AC_PREFIX_DEFAULT([`pwd`])

#############################################################################
#                         Standard build tool stuff                         #
#############################################################################

# Get the system type
AC_CANONICAL_BUILD

# If this project depends on external projects, the Externals file in
# the source root directory contains definition of where to find those
# externals.  The following macro ensures that those externals are
# retrieved by svn if they are not there yet.
AC_COIN_PROJECTDIR_INIT(Smi)

# Check if user wants to produce debugging code
AC_COIN_DEBUG_COMPILE(Smi)

# Get the name of the C++ compiler and appropriate compiler options
AC_COIN_PROG_CXX

# Initialize automake and libtool
AC_COIN_INIT_AUTO_TOOLS

#############################################################################
#                              COIN components                              #
#############################################################################

AC_COIN_CHECK_PACKAGE(Osi, [osi], [Smi])
if test $coin_has_osi != yes ; then
  AC_MSG_ERROR([Required package Osi not found.])
fi
AC_COIN_CHECK_PACKAGE(Clp, [osi-clp])
AC_COIN_CHECK_PACKAGE(DataStochastic, [coindatastochastic])
AC_COIN_CHECK_PACKAGE(FlopCpp, [flopcpp])

##############################################################################
#                   VPATH links for example input files                      #
##############################################################################

AC_COIN_VPATH_LINK(test/SmiTestData/bug.cor)
AC_COIN_VPATH_LINK(test/SmiTestData/bug.stoch)
AC_COIN_VPATH_LINK(test/SmiTestData/bug.time)
AC_COIN_VPATH_LINK(test/SmiTestData/app0110.cor)
AC_COIN_VPATH_LINK(test/SmiTestData/app0110.stoch)
AC_COIN_VPATH_LINK(test/SmiTestData/app0110.time)
AC_COIN_VPATH_LINK(test/SmiTestData/app0110R.cor)
AC_COIN_VPATH_LINK(test/SmiTestData/app0110R.stoch)
AC_COIN_VPATH_LINK(test/SmiTestData/app0110R.time)

# In case this is a VPATH configuration we need to make sure that the
# input files for the examples are available in the VPATH directory.

##############################################################################
#                   Finishing up by writing all the output                   #
##############################################################################

# Here list all the files that configure should create (except for the
# configuration header file)
AC_CONFIG_FILES([Makefile
                 examples/Makefile
                 src/Makefile
                 test/Makefile
                 flopcpp_examples/Makefile
                 smi.pc
                 smi-uninstalled.pc])

# Here put the location and name of the configuration header file
# Smi does not make use of the config header files, so the first is
# just a dummy file. The second is for a user who wants to get the
# Smi version number.
AC_CONFIG_HEADER([src/config.h src/config_smi.h])

# Finally, we let configure write all the output...
AC_COIN_FINALIZE
