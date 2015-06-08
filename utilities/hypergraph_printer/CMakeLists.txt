file(GLOB_RECURSE PRINTER_SRCS *.cpp *.cxx *.c)
list(APPEND PRINTER_SRCS "../shared/src/StringUtils.cpp")

include_directories("include")
include_directories("../shared/include")

set(PRINTER_BIN ${PROJECT_NAME}-hypergraph-printer)
add_executable(${PRINTER_BIN} ${PRINTER_SRCS})

# Require that it is compiled with C++11.
set_property(TARGET ${PRINTER_BIN} PROPERTY CXX_STANDARD 11)
set_property(TARGET ${PRINTER_BIN} PROPERTY CXX_STANDARD_REQUIRED ON)
