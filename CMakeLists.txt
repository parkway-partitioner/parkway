# Minimum CMake version.
cmake_minimum_required(VERSION 3.2)


### Project information #######################################################
# Application information and version.
project("parkway")
set(APPLICATION_NAME "Parkway")
set(VERSION_MAJOR 2)
set(VERSION_MINOR 12)
set(VERSION_STRING "${VERSION_MAJOR}.${VERSION_MINOR}")


### Output directories ########################################################
set(MAINFOLDER ${PROJECT_SOURCE_DIR})
set(EXECUTABLE_OUTPUT_PATH "${MAINFOLDER}/bin")
set(LIBRARY_OUTPUT_PATH "${MAINFOLDER}/lib")

# Set the name of the output library.
set(PROJECT_LIB ${PROJECT_NAME})

set_directory_properties(PROPERTIES EP_PREFIX ${CMAKE_BINARY_DIR}/thirdparty)

### Extra CMake opptions ######################################################
# Add a CMake build type requirement. Default to release with debug info.
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "RelWithDebInfo" CACHE STRING
      "Choose the type of build, options are: Debug Developer MinSizeRel
      Release RelWithDebInfo." FORCE)
endif()

if(NOT PARKWAY_TESTS)
  set(PARKWAY_TESTS false CACHE BOOL
      "Choose to build with tests, options are: true false." FORCE)
endif()

# Add -Wall and -pedantic if GCC / Clang and building in debug mode.
if(${CMAKE_BUILD_TYPE} STREQUAL "debug")
  if(${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU" OR
     ${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang")
    message(STATUS "Configuring in debug mode, adding '-Wall -pedantic' flags")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -pedantic")
  endif()
endif()


### Packages ##################################################################
# Find any required packages.
find_package(MPI REQUIRED)
find_package(Boost COMPONENTS log_setup log thread system REQUIRED)
add_definitions(-DBOOST_LOG_DYN_LINK)



### Testing framework #########################################################
# Check if enabled.
if(PARKWAY_TESTS)
  # Enable CTest to handle running tests with 'make test'.
  include(CTest)

  add_subdirectory("${MAINFOLDER}/tests/unit_tests")
endif()


### Include paths #############################################################
# parkway
include_directories("${MAINFOLDER}/include")
# MPI
include_directories(${MPI_INCLUDE_PATH})
# Boost
include_directories(${Boost_INCLUDE_PATH})


### CMakeLists to consider ####################################################
# libparkway
add_subdirectory("${MAINFOLDER}/src")

# utilities
add_subdirectory("${MAINFOLDER}/utilities/hypergraph_converter")
add_subdirectory("${MAINFOLDER}/utilities/hypergraph_printer")
add_subdirectory("${MAINFOLDER}/utilities/driver")
