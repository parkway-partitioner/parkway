file(GLOB_RECURSE PRINTER_SRCS *.cpp *.cxx *.c)
file(GLOB_RECURSE PRINTER_HEADERS include/*.h include/*.hpp
  ../shared/include/*.h ../shared/include/*.hpp)
list(APPEND PRINTER_SRCS "../shared/src/StringUtils.cpp")

include_directories("include")
include_directories("../shared/include")

set(PRINTER_BIN ${PROJECT_NAME}-hypergraph-printer)
# Include headers here so that generated IDE projects include them.
add_executable(${PRINTER_BIN} ${PRINTER_SRCS} ${PRINTER_HEADERS})

# Require that it is compiled with C++11.
set_property(TARGET ${PRINTER_BIN} PROPERTY CXX_STANDARD 11)
set_property(TARGET ${PRINTER_BIN} PROPERTY CXX_STANDARD_REQUIRED ON)
