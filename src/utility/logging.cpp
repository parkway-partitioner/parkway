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

namespace status {

stream::stream(const char *filename)
  : filestream_(filename), stream_(filestream_), buffer_(&null_buffer_) {
}

stream::stream() : stream_(std::cout), buffer_(&null_buffer_) {
}

void stream::enable() {
  enabled_ = true;
}

void stream::disable() {
  enabled_ = false;
}

std::ostream &stream::operator()(const level_t lvl) {
  if (enabled_ && lvl <= filter_) {
    return stream_;
  }
  return buffer_;
}

void stream::set_filter_level(const level_t lvl) {
  filter_ = lvl;
}


status_handler::status_handler() : output_stream_(new stream) {
}


void status_handler::disable() {
  output_stream_->disable();
}

void status_handler::enable() {
  output_stream_->enable();
}

void status_handler::set_rank(std::size_t rank) {
  rank_ = rank;
  if (rank_ == valid_rank_) {
    enable();
  } else {
    disable();
  }
}

std::ostream &status_handler::operator()(const level_t lvl) {
  return static_cast<std::ofstream&>(output_stream_->operator()(lvl));
}


void status_handler::set_filter_level(const level_t lvl) {
  output_stream_->set_filter_level(lvl);
}

void status_handler::set_output_file(const char *filename) {
  output_stream_.release();
  output_stream_.reset(new stream(filename));
}

bool stream::enabled_ = true;
level_t stream::filter_ = level_t::progress;
std::size_t status_handler::valid_rank_ = 0;
std::size_t status_handler::rank_ = 0;

}
}
}
