#include "utility/logging.hpp"

#include "boost/log/expressions.hpp"
#include "boost/log/sinks/text_file_backend.hpp"
#include "boost/log/sources/record_ostream.hpp"
#include "boost/log/utility/setup/common_attributes.hpp"
#include "boost/log/utility/setup/file.hpp"

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

void set_log_file(const std::string &basename, const std::size_t rank) {
  // Add attributes used by the formatter (e.g. TimeStamp, etc)
  boost::log::add_common_attributes();

  // Add Severity level (e.g. 'info', 'error', etc.) to the formatter
  boost::log::register_simple_formatter_factory
    <logging::trivial::severity_level, char>("Severity");

  std::string rank_str = boost::lexical_cast<std::string>(rank);

  boost::log::add_file_log(
    keywords::file_name = basename + "_p" + rank_str + "__%Y-%m-%d__%H%M.log",
    keywords::format = "[" + rank_str + "][%LineID%][%Severity%]: %Message%");
}

}  // namespace logging
}  // namespace utility
}  // namespace parkway

namespace parkway {
namespace utility {
namespace status {

std::ostream *handler::out_ = &std::cout;
std::ofstream *handler::filestream_ = nullptr;
std::unique_ptr<char> handler::buffer_(new char[1024]);

std::size_t handler::write_on_ = 0;
std::size_t handler::rank_ = 0;
bool handler::show_info_ = true;
bool handler::show_progress_ = true;

// This status handler exists to call tidy up when the program exits.
handler __handler__;

}  // namespace status
}  // namespace utility
}  // namespace parkway
