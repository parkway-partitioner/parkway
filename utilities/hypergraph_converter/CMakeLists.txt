file(GLOB_RECURSE CONVERTER_SRCS *.cpp *.cxx *.c)
list(APPEND CONVERTER_SRCS "../shared/src/StringUtils.cpp")

include_directories("include")
include_directories("../shared/include")

set(CONVERTER_BIN ${PROJECT_NAME}-hypergraph-converter)
add_executable(${CONVERTER_BIN} ${CONVERTER_SRCS})

# Require that it is compiled with C++11.
set_property(TARGET ${CONVERTER_BIN} PROPERTY CXX_STANDARD 11)
set_property(TARGET ${CONVERTER_BIN} PROPERTY CXX_STANDARD_REQUIRED ON)
