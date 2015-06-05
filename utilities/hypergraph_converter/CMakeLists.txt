file(GLOB_RECURSE CONVERTER_SRCS *.cpp *.cxx *.c)
list(APPEND CONVERTER_SRCS "../shared/src/StringUtils.cpp")

include_directories("include")
include_directories("../shared/include")

set(CONVERTER_BIN ${PROJECT_NAME}-hypergraph-converter)
add_executable(${CONVERTER_BIN} ${CONVERTER_SRCS})
