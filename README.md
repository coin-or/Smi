# Smi - Stochastic Modeling Interface

To install this package, please look at the INSTALL file.

This package contains several subdirectories corresponding to COIN-OR
projects (www.coin-or.org). The AUTHORS, LICENSE and README files in 
each of the subdirectories give more information about these projects.

The examples directory contains a few use cases.

* stoch.cpp is the main repository of use cases -- it contains several
  distinct examples of input styles.

* investment.cpp is a use case using FlopC++.  This is under active development
  but is stable enough to be included in trunk

## SMI FAQs

### What is SMI?
SMI stands for Stochastic Modeling Interface. It is an interface for problems in which uncertainty and optimization appear together. There are many modeling and algorithmic approaches that could belong here, like: recourse programming, chance constrained programming, stochastic control and dynamic programming, robust optimization, etc, etc. SMI is intended to be like OSI in the sense that an SmiXX object is an implementation derived from a base class that takes care of a number of commonly encountered programming issues, like handling probability distributions, managing problem generation, interacting with solvers to obtain solution information, etc.


### What is in the current release?
The current release implements a multiperiod scenario stochastic programming object called SmiScnModel. It supports an SMPS file reader method, a direct "genScenario" method, a method to generate a deterministic equivalent, and several methods to get solution data by scenario. This is a fully native COIN-OR implementation.


### Do I need to install any particular solver?
No. It can use any OSI-compatible solver.


### How do I use it?
You can obtain the Smi source code either via subversion or in form of nightly generated tarballs.  The recommended method is to use subversion because it makes it easier to obtain updates. The following commands may be used to obtain and build Smi from the source code using subversion:

 1. `svn co https://projects.coin-or.org/svn/Smi/stable/0.96 smi_stable`
 1. `cd smi_stable`
 1. `./configure -C`
 1. `make`
 1. `make test`
 1. `make install`
 1. `cd Smi/examples`
 1. `make`

The last step compiles the file Smi/example/stoch.cpp which has a number of use cases. 

If the "make test" fails, it is likely because the SMPS files were not found.  In the
top of the unitTest.cpp file there is a line:
`#define SMI_TEST_DATA_DIR "./SmiTestData"`

The make test step will fail to find the SMPS files if this directory does not resolve to `<download_root>/coin-Smi_0.96/Smi/test/SmiTestData`.

WINDOWS users: The [BuildTools](http://projects.coin-or.org/BuildTools/wiki) project has additional details on downloading, building, installing (e.g. for Windows), available options and troubleshooting.  Note -- you may have to open the Microsoft Visual Studio properties for projects `examples` and `unitTestSmi` and check to see that the working directory is set to:  `$(SolutionDir)/../../test`


## Included Projects

If you download the Smi package, you get [these](Dependencies) additional projects.
