file(GLOB_RECURSE DRIVER_SRCS *.cpp *.cxx *.c)
set(PROJECT_LIBRARIES ${MPI_LIBRARIES} ${PROJECT_LIB})
include_directories("include")

set(DRIVER_BIN ${PROJECT_NAME}-driver)
add_executable(${DRIVER_BIN} ${DRIVER_SRCS})
target_link_libraries(${DRIVER_BIN} ${PROJECT_LIBRARIES})
