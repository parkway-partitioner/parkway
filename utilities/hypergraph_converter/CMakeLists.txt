file(GLOB_RECURSE CONVERTER_SRCS *.cpp *.cxx *.c)
file(GLOB_RECURSE CONVERTER_HEADERS include/*.h include/*.hpp
  ../shared/include/*.h ../shared/include/*.hpp)
list(APPEND CONVERTER_SRCS "../shared/src/StringUtils.cpp")

include_directories("include")
include_directories("../shared/include")

set(CONVERTER_BIN ${PROJECT_NAME}-hypergraph-converter)
# Include headers here so that generated IDE projects include them.
add_executable(${CONVERTER_BIN} ${CONVERTER_SRCS} ${CONVERTER_HEADERS})

# Require that it is compiled with C++11.
set_property(TARGET ${CONVERTER_BIN} PROPERTY CXX_STANDARD 11)
set_property(TARGET ${CONVERTER_BIN} PROPERTY CXX_STANDARD_REQUIRED ON)
