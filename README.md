# Par<i>k</i>way

Parallel hypergraph partitioning in C++ using MPI.

## Minimum Requirements

* [CMake](http://cmake.org/) 3.2
* C++11 compliant compiler (tested with Clang 3.6, GCC 4.8.2, GCC 5.1.0)
* MPI (tested with MPICH 3.0.4 and MPICH 3.1.4)

## Building

Run the `configure` script with any additional CMake arguments (i.e. to select
a generator). Parkway options are as such:

* `-D PARKWAY_TESTS=<true|false>` to add support for building and running unit
  tests.

Assuming a UNIX system using GNU Makefiles, once configured `make help` will
list all build targets. A number of targets are:

* `make parkway` - builds `libparkway`.
* `make parkway_driver` - builds a utility to partition hypergraphs using
  `libparkway`.
* `make` - builds the above and a few more utilities.
