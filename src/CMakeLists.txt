# Find source files.
file(GLOB_RECURSE PROJECT_SRCS *.cpp *.cxx *.c)
file(GLOB_RECURSE PROJECT_HEADERS
  "${MAINFOLDER}/include/*.h" "${MAINFOLDER}/include/*.hpp")
set(PROJECT_LIBRARIES ${MPI_LIBRARIES})

# Compile the sources into a library and link with 3rd party libraries.
# Include headers here so that generated IDE projects include them.
add_library(${PROJECT_LIB} ${PROJECT_SRCS} ${PROJECT_HEADERS})

# Require that libparkway is compiled with C++11.
set_property(TARGET ${PROJECT_LIB} PROPERTY CXX_STANDARD 11)
set_property(TARGET ${PROJECT_LIB} PROPERTY CXX_STANDARD_REQUIRED ON)

# Link parkway to the required libraries.
target_link_libraries(${PROJECT_LIB} ${PROJECT_LIBRARIES})
