# Minimum CMake version.
cmake_minimum_required(VERSION 2.6)

# Application information and version.
project("parkway")
set(APPLICATION_NAME "parkway")
set(VERSION_MAJOR 2)
set(VERSION_MINOR 12)
set(VERSION_STRING "${VERSION_MAJOR}.${VERSION_MINOR}")


# Output directories.
set(MAINFOLDER ${PROJECT_SOURCE_DIR})
set(EXECUTABLE_OUTPUT_PATH "${MAINFOLDER}/bin")
set(LIBRARY_OUTPUT_PATH "${MAINFOLDER}/lib")

# Output binary and library names to "parkway" and "libparkway".
set(PROJECT_BIN ${PROJECT_NAME})
set(PROJECT_LIB "lib${PROJECT_NAME}")

# Find any required packages.
find_package(MPI REQUIRED)

# Look for headers in "include".
include_directories("${MAINFOLDER}/include")
# Look for MPI headers.
include_directories(${MPI_INCLUDE_PATH})

# Look for sources.
add_subdirectory("${MAINFOLDER}/src")

# Look for utilities.
add_subdirectory("${MAINFOLDER}/utilities/hypergraph_converter")
add_subdirectory("${MAINFOLDER}/utilities/hypergraph_printer")
add_subdirectory("${MAINFOLDER}/utilities/driver")
