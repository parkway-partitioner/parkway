#ifndef UTILITY_LOGGING_HPP_
#define UTILITY_LOGGING_HPP_
#include <string>
#include "boost/log/trivial.hpp"

#define LOG BOOST_LOG_TRIVIAL

namespace parkway {
namespace utility {
namespace logging {

// Sets the maximum reported severity level.
// Options are "all", "trace", "debug", "info", "warn", "error", "fatal", and
// "off".
void set_filter_level(const std::string &level);

void set_log_file(const std::string &filename, const std::size_t rank);

}
}
}


#endif  // UTILITY_LOGGING_HPP_
