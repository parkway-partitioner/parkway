file(GLOB_RECURSE DRIVER_SRCS *.cpp *.cxx *.c)
set(PROJECT_LIBRARIES ${MPI_LIBRARIES} ${PROJECT_LIB})
include_directories("include")

set(DRIVER_BIN ${PROJECT_NAME}-driver)
add_executable(${DRIVER_BIN} ${DRIVER_SRCS})

# Require that it is compiled with C++11.
set_property(TARGET ${DRIVER_BIN} PROPERTY CXX_STANDARD 11)
set_property(TARGET ${DRIVER_BIN} PROPERTY CXX_STANDARD_REQUIRED ON)

target_link_libraries(${DRIVER_BIN} ${PROJECT_LIBRARIES})
