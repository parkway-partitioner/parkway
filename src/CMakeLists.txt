# Find source files.
file(GLOB_RECURSE PROJECT_SRCS *.cpp *.cxx *.c)
set(PROJECT_LIBRARIES ${MPI_LIBRARIES} ${Boost_LIBRARIES})

# Compile the sources into a library and link with 3rd party libraries.
add_library(${PROJECT_LIB} ${PROJECT_SRCS})

# Require that libparkway is compiled with C++11.
set_property(TARGET ${PROJECT_LIB} PROPERTY CXX_STANDARD 11)
set_property(TARGET ${PROJECT_LIB} PROPERTY CXX_STANDARD_REQUIRED ON)

# Link parkway to the required libraries.
target_link_libraries(${PROJECT_LIB} ${PROJECT_LIBRARIES})
