#include "utility/logging.hpp"
#include <boost/log/expressions.hpp>
#include <boost/log/sources/record_ostream.hpp>
#include <boost/log/utility/setup/common_attributes.hpp>

namespace parkway {
namespace utility {
namespace logging {

namespace trivial = boost::log::trivial;
namespace keywords = boost::log::keywords;

void set_filter_level(std::string const &level) {
  // The return type depends on whether compilation is whith threads.
  auto logger = boost::log::core::get();

  if (level == "all" || level == "trace") {
    logger->set_filter(trivial::severity >= trivial::trace);
    LOG(info) << "Log level set to TRACE";
  } else if (level == "debug") {
    logger->set_filter(trivial::severity >= trivial::debug);
    LOG(info) << "Log level set to DEBUG";
  } else if (level == "info") {
    logger->set_filter(trivial::severity >= trivial::info);
    LOG(info) << "Log level set to INFO";
  } else if (level == "warn") {
    logger->set_filter(trivial::severity >= trivial::warning);
  } else if (level == "error") {
    logger->set_filter(trivial::severity >= trivial::error);
  } else if (level == "fatal") {
    logger->set_filter(trivial::severity >= trivial::fatal);
  } else if ((level == "off") || (level == "none")) {
    logger->set_logging_enabled(false);
  } else {
    logger->set_filter(trivial::severity >= trivial::fatal);
  }
}

}
}
}
