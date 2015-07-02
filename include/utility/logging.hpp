#ifndef UTILITY_LOGGING_HPP_
#define UTILITY_LOGGING_HPP_
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include "boost/log/trivial.hpp"
#include "utility/status.hpp"

#define LOG BOOST_LOG_TRIVIAL

using parkway::utility::status::info;
using parkway::utility::status::progress;

namespace parkway {
namespace utility {
namespace logging {

// Sets the maximum reported severity level.
// Options are "all", "trace", "debug", "info", "warn", "error", "fatal", and
// "off".
void set_filter_level(const std::string &level);

void set_log_file(const std::string &filename, const std::size_t rank);

}  // namespace logging
}  // namespace utility
}  // namespace parkway


#endif  // UTILITY_LOGGING_HPP_
