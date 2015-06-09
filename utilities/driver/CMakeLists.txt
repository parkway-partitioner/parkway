file(GLOB_RECURSE DRIVER_SRCS *.cpp *.cxx *.c)
file(GLOB_RECURSE DRIVER_HEADERS include/*.h include/*.hpp)
set(PROJECT_LIBRARIES ${MPI_LIBRARIES} ${PROJECT_LIB})
include_directories("include")

set(DRIVER_BIN ${PROJECT_NAME}_driver)
# Include headers here so that generated IDE projects include them.
add_executable(${DRIVER_BIN} ${DRIVER_SRCS} ${DRIVER_HEADERS})

# Require that it is compiled with C++11.
set_property(TARGET ${DRIVER_BIN} PROPERTY CXX_STANDARD 11)
set_property(TARGET ${DRIVER_BIN} PROPERTY CXX_STANDARD_REQUIRED ON)

target_link_libraries(${DRIVER_BIN} ${PROJECT_LIBRARIES})
