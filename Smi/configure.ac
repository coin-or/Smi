# Copyright (C) 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.

#############################################################################
#                       Names and other basic things                        #
#############################################################################

AC_INIT([Smi],[devel],[https://github.com/coin-or/Smi/issues/new],[coin-or-smi],[https://github.com/coin-or/Smi])

AC_COPYRIGHT([
Copyright 2006 International Business Machines and others.
All Rights Reserved.
This file is part of the open source package Smi which is distributed
under the Eclipse Public License.])

# List one file in the package so that the configure script can test
# whether the package is actually there
AC_CONFIG_SRCDIR(src/SmiScnModel.hpp)

# Do some initialization work (version numbers, change prefix default, ...)
AC_COIN_INITIALIZE

#############################################################################
#                         Standard build tool stuff                         #
#############################################################################

# Get the name of the C++ compiler and appropriate compiler options
AC_COIN_PROG_CXX

# This is a C++ package, set the language accordingly.
AC_LANG_PUSH(C++)

# Initialize libtool
AC_COIN_PROG_LIBTOOL

# set RPATH_FLAGS to the compiler link flags required to hardcode location
# of the shared objects (expanded_libdir is set somewhere in configure before)
# (used in the examples)
AC_COIN_RPATH_FLAGS([$expanded_libdir])

#############################################################################
#                              COIN-OR components                           #
#############################################################################

AC_COIN_CHK_PKG(Osi,[SmiLib SmiTest])
if test $coin_has_osi != yes ; then
  AC_MSG_ERROR([Required package Osi not found.])
fi

# The test and the examples use Clp (one .pc file per check, as some
# pkg-config versions do not accept several modules in one argument)
AC_COIN_CHK_PKG(Clp,[SmiTest])
AC_COIN_CHK_PKG(OsiClp,[SmiTest],[osi-clp])

# FlopC++ is used by the examples in flopcpp_examples
AC_COIN_CHK_PKG(FlopCpp,,[flopcpp])

# Stochastic programming instances for the examples
AC_COIN_CHK_PKG(DataStochastic,,[coindatastochastic],,dataonly)

##############################################################################
#                   VPATH links for test input files                         #
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

##############################################################################
#                   Finishing up by writing all the output                   #
##############################################################################

# Doxygen documentation
AC_COIN_DOXYGEN([CoinUtils Osi Clp Cgl FlopCpp])

# Finalize libs and cflags.
AC_COIN_FINALIZE_FLAGS([SmiLib SmiTest])

# Here list all the files that configure should create (except for the
# configuration header file)
AC_CONFIG_FILES([Makefile
                 examples/Makefile
                 src/Makefile
                 test/Makefile
                 flopcpp_examples/Makefile
                 smi.pc])

AC_CONFIG_FILES([doxydoc/doxygen.conf])

# Here put the location and name of the configuration header file
AC_CONFIG_HEADERS([src/config.h src/config_smi.h])

# Finally, we let configure write all the output...
AC_COIN_FINALIZE
