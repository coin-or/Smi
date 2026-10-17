# Copyright (C) 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.

include $(top_srcdir)/BuildTools/Makemain.inc

if COIN_HAS_CLP

noinst_PROGRAMS = unitTest

unitTest_SOURCES = unitTest.cpp

unitTest_LDADD = ../src/libSmi.la $(SMITEST_LFLAGS)

AM_CPPFLAGS = -I$(srcdir)/../src -I../src $(SMITEST_CFLAGS)
AM_LDFLAGS = $(LT_LDFLAGS)

test: unitTest$(EXEEXT)
	./unitTest$(EXEEXT)

else

test:
	echo "Clp required to run unittest."

endif

.PHONY: test

DISTCLEANFILES = bug_gen.mps
