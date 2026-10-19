# Smi 0.96

[![A COIN-OR Project](https://coin-or.github.io/coin-or-badge.png)](https://www.coin-or.org)

Projects such as this one are maintained by a small group of volunteers under
the auspices of the non-profit [COIN-OR Foundation](https://www.coin-or.org)
and we need your help! Please consider [sponsoring our
activities](https://github.com/sponsors/coin-or) or [volunteering](mailto:volunteer@coin-or.org) to help!

[![Latest Release](https://img.shields.io/github/v/release/coin-or/Smi?sort=semver)](https://github.com/coin-or/Smi/releases)

_This file is auto-generated from [config.yml](.coin-or/config.yml) using the 
[generate_readme](.coin-or/generate_readme) script.
To make changes, please edit [config.yml](.coin-or/config.yml) or the generation scripts
[here](.coin-or/generate_readme) and [here](https://github.com/coin-or/coinbrew/blob/master/scripts/generate_readme)._

Smi (Stochastic Modeling Interface) is an interface for problems in which
uncertainty and optimization appear together. It is intended to be like
[Osi](https://github.com/coin-or/Osi) in the sense that an `SmiXX` object is
derived from a base class that takes care of a number of commonly encountered
programming issues, like handling probability distributions, managing problem
generation, and interacting with solvers to obtain solution information.

Smi implements a multiperiod scenario stochastic programming object called
`SmiScnModel`. It supports
 * a reader for Stochastic MPS (SMPS) files,
 * generating scenario trees directly from paths and from discrete random
   variables,
 * generating the deterministic equivalent problem, which can be solved by
   any Osi-compatible solver, and
 * getting solution data by stage and scenario.

Smi is written in C++ and is released as open source under the [Eclipse Public License 2.0](http://www.opensource.org/licenses/EPL-2.0).

It is distributed under the auspices of the [COIN-OR Foundation](https://www.coin-or.org).

The Smi development site is https://github.com/coin-or/Smi.

## CURRENT BUILD STATUS

[![Windows Builds](https://github.com/coin-or/Smi/actions/workflows/windows-ci.yml/badge.svg?branch=stable/0.96)](https://github.com/coin-or/Smi/actions/workflows/windows-ci.yml?query=branch%3Astable/0.96)

[![Linux and MacOS Builds](https://github.com/coin-or/Smi/actions/workflows/linux-ci.yml/badge.svg?branch=stable/0.96)](https://github.com/coin-or/Smi/actions/workflows/linux-ci.yml?query=branch%3Astable/0.96)

## DOWNLOAD

### Docker image

There is a Docker image that provides Smi, as well as other projects
in the [COIN-OR Optimization
Suite](https://github.com/coin-or/COIN-OR-OptimizationSuite) [here](https://hub.docker.com/r/coinor/coin-or-optimization-suite).

### Binaries

For newer releases, binaries will be made available as assets attached to
releases in Github
[here](https://github.com/coin-or/Smi/releases). Older binaries
are archived as part of Smi
[here](https://www.coin-or.org/download/binary/Smi).

Due to license incompatibilities, pre-compiled binaries may lack some
functionality. If binaries are not available for your platform for the latest
version and you would like to request them to be built and posted, feel free
to open a discussion on Github.

### Source

Source code can be obtained either by

 * Downloading a snapshot of the source code for the latest release version of Smi from the
 [releases](https://github.com/coin-or/Smi/releases) page.
 * Cloning this repository from [Github](https://github.com/coin-or/Smi) or 
 * Using the [coinbrew](https://github.com/coin-or/coinbrew) script to get the project and all dependencies (recommended, see below).   

Below is a quick start guide for building on common platforms. More detailed
build instructions are
[here](https://coin-or.github.io/user_introduction.html).

### Dependencies

Smi has a number of dependencies, which are detailed in
[config.yml](.coin-or/config.yml). Dependencies on other COIN-OR projects are
automatically downloaded when obtaining the source with `coinbrew`. For some
of the remaining third-party dependencies, automatic download scripts and
build wrappers are provided (and will also be automatically run for required
and recommended dependencies), while other libraries that are easy to obtain
must be installed using an appropriate package manager (or may come with your
OS by default). 

## BUILDING from source

The quick start assumes you are in a bash shell. 

### Using `coinbrew`

To download and build Smi from source, execute the 
following on the command line. 
```
wget https://raw.githubusercontent.com/coin-or/coinbrew/master/coinbrew
chmod u+x coinbrew
./coinbrew fetch Smi@0.96
./coinbrew build Smi
```
For more detailed instructions on coinbrew, see https://coin-or.github.io/coinbrew.
The `coinbrew` script will fetch the additional projects specified in the Dependencies section of [config.yml](.coin-or/config.yml).

### Without `coinbrew` (Expert users)

 * Download the source code, e.g., by cloning the git repo https://github.com/coin-or/Smi
 * Download and install the source code for the dependencies listed in [config.yml](.coin-or/config.yml)
 * Build the code as follows (make sure to set `PKG_CONFIG_PATH` to the `lib/pkgconfig` directory of the installation of the dependencies).

```
./configure -C
make
make test
make install
```

### Examples

The `examples` directory has a number of use cases. `stoch.cpp` contains
several distinct examples of input styles; it reads instances from
[Data-Stochastic](https://github.com/coin-or-tools/Data-Stochastic) if that is
available. The `flopcpp_examples` directory has examples that use
[FlopC++](https://github.com/coin-or/FlopCpp). After building and installing
Smi, they can be built by running `make` in the corresponding subdirectory of
the build directory, e.g.,
```
cd build/Smi/0.96/examples
make
./stoch
```

## Doxygen Documentation

If you have `Doxygen` available, you can build a HTML documentation by typing

`make doxygen-docs` 

in the build directory. If Smi was built via `coinbrew`, then the build
directory will be `./build/Smi/0.96` by default. The doxygen documentation main file
is found at `<build-dir>/doxydoc/html/index.html`.


## Project Links

 * [Code of Conduct](https://www.coin-or.org/code-of-conduct/)
 * [COIN-OR Web Site](http://www.coin-or.org/)
 * [COIN-OR general discussion forum](https://github.com/orgs/coin-or/discussions)
 * [Smi Discussion forum](https://github.com/coin-or/Smi/discussions)
 * [Report a bug](https://github.com/coin-or/Smi/issues/new)

