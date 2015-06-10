# GMock (and GTest) require threads.
find_package(Threads REQUIRED)

# GMock (and GTest) will be downloaded as external project.
include(ExternalProject)

# Add GMock as an external project (this will download the latest version
# from the SVN repo and build it in 'thirdparty'.)
ExternalProject_Add(
  googlemock
  SVN_REPOSITORY http://googlemock.googlecode.com/svn/trunk/
  TIMEOUT 30
  # Force separate output paths for debug and release builds to allow easy
  # identification of correct lib in subsequent TARGET_LINK_LIBRARIES
  # commands.
  # Disable install step.
  INSTALL_COMMAND ""
  # Wrap download, configure and build steps in a script to log output.
  LOG_DOWNLOAD ON
  LOG_CONFIGURE ON
  LOG_BUILD ON
)

# Find the test source files.
file(GLOB_RECURSE UNIT_TEST_FILES *.cpp *.cxx *.c)

# Set the name of the output test driver and the dependencies.
set(PROJECT_TEST_DRIVER "${PROJECT_NAME}_test_driver")
add_executable(${PROJECT_TEST_DRIVER} ${UNIT_TEST_FILES})
add_dependencies(${PROJECT_TEST_DRIVER} googlemock)
include_directories("${MAINFOLDER}/include")


# Require that libparkway is compiled with C++11.
set_property(TARGET ${PROJECT_TEST_DRIVER} PROPERTY CXX_STANDARD 11)
set_property(TARGET ${PROJECT_TEST_DRIVER} PROPERTY CXX_STANDARD_REQUIRED ON)

# Specify include directories for googlemock and googletest.
ExternalProject_Get_Property(googlemock source_dir binary_dir)
include_directories("${source_dir}/include")
include_directories("${source_dir}/gtest/include")
target_link_libraries(${PROJECT_TEST_DRIVER} ${PROJECT_LIB}
  ${binary_dir}/${CMAKE_FIND_LIBRARY_PREFIXES}gmock.a
  ${binary_dir}/${CMAKE_FIND_LIBRARY_PREFIXES}gmock_main.a
  ${binary_dir}/gtest/${CMAKE_FIND_LIBRARY_PREFIXES}gtest.a
  ${binary_dir}/gtest/${CMAKE_FIND_LIBRARY_PREFIXES}gtest_main.a
  ${CMAKE_THREAD_LIBS_INIT}
)

# Scan the test sources to find test name and add them to the CTest.
foreach(UNIT_TEST_FILE ${UNIT_TEST_FILES})
  message(STATUS "Test file: ${UNIT_TEST_FILE} ")
  file(STRINGS ${UNIT_TEST_FILE} GTEST_NAMES REGEX ^TEST)
  foreach(GTEST_NAME ${GTEST_NAMES})
    string(REGEX REPLACE ["\) \(,"] ";" GTEST_NAME ${GTEST_NAME})
    list(GET GTEST_NAME 1 GTEST_GROUP_NAME)
    list(GET GTEST_NAME 3 GTEST_NAME)
    add_test(
      NAME ${GTEST_GROUP_NAME}.${GTEST_NAME}
      COMMAND "${EXECUTABLE_OUTPUT_PATH}/${PROJECT_TEST_DRIVER}"
              "--gtest_filter=${GTEST_GROUP_NAME}.${GTEST_NAME}"
    )
  endforeach()
endforeach()

# Add a target to check that tests have been built.
add_custom_target(check
  COMMAND ${CMAKE_CTEST_COMMAND}
  DEPENDS ${PROJECT_TEST_DRIVER}
)
