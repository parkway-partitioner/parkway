#ifndef UTILITY_LOGGING_HPP_
#define UTILITY_LOGGING_HPP_
#include "boost/log/trivial.hpp"

#define LOG BOOST_LOG_TRIVIAL

namespace parkway {
namespace utility {
namespace logging {

// Sets the maximum reported severity level.
// Options are "all", "trace", "debug", "info", "warn", "error", "fatal", and
// "off".
void set_filter_level(std::string const &level);


}
}
}


#endif  // UTILITY_LOGGING_HPP_
