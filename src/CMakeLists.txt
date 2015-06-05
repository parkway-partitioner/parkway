# Find source files.
file(GLOB_RECURSE PROJECT_SRCS *.cpp *.cxx *.c)
set(PROJECT_LIBRARIES ${MPI_LIBRARIES})

# Compile the sources into a library and link with 3rd party libraries.
add_library(${PROJECT_LIB} ${PROJECT_SRCS})
target_link_libraries(${PROJECT_LIB} ${PROJECT_LIBRARIES})
