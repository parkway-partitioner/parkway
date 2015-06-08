# Minimum CMake version.
cmake_minimum_required(VERSION 2.6)


# Application information and version.
project("parkway")
set(APPLICATION_NAME "Parkway")
set(VERSION_MAJOR 2)
set(VERSION_MINOR 12)
set(VERSION_STRING "${VERSION_MAJOR}.${VERSION_MINOR}")

# Output directories.
set(MAINFOLDER ${PROJECT_SOURCE_DIR})
set(EXECUTABLE_OUTPUT_PATH "${MAINFOLDER}/bin")
set(LIBRARY_OUTPUT_PATH "${MAINFOLDER}/lib")

# Output binary and library names to "parkway" and "libparkway".
set(PROJECT_BIN ${PROJECT_NAME})
set(PROJECT_LIB ${PROJECT_NAME})

# Add a CMake build type requirement. Default to release with debug info.
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "RelWithDebInfo" CACHE STRING
      "Choose the type of build, options are: Debug Developer MinSizeRel
      Release RelWithDebInfo." FORCE)
endif()

# Add -Wall and -pedantic if GCC / Clang and building in debug mode.
if(${CMAKE_BUILD_TYPE} STREQUAL "debug")
  if(${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU" OR
     ${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang")
    message(STATUS "Configuring in debug mode, adding '-Wall -pedantic' flags")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -pedantic")
  endif()
endif()

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
