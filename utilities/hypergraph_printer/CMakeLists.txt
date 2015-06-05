file(GLOB_RECURSE PRINTER_SRCS *.cpp *.cxx *.c)
list(APPEND PRINTER_SRCS "../shared/src/StringUtils.cpp")

include_directories("include")
include_directories("../shared/include")

set(PRINTER_BIN ${PROJECT_NAME}-hypergraph-printer)
add_executable(${PRINTER_BIN} ${PRINTER_SRCS})
