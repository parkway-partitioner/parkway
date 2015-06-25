#ifndef UTILITY_LOGGING_HPP_
#define UTILITY_LOGGING_HPP_
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
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

namespace status {

enum class level_t {info = 0, progress};

class null_buffer : public std::streambuf {
 public:
  int overflow(int c) {
    return c;
  }
};

class stream {
 public:
  stream();
  stream(const char *filename);

  void disable();
  void enable();

  std::ostream &operator()(const level_t lvl);
  void set_filter_level(const level_t lvl);


 private:
  static bool enabled_;
  std::ofstream filestream_;
  std::ostream &stream_;
  null_buffer null_buffer_;
  std::ostream buffer_;
  static level_t filter_;
};

class status_handler {
 public:
  status_handler();

  void disable();
  void enable();

  void set_rank(std::size_t rank);
  void set_filter_level(const level_t lvl);
  void set_output_file(const char *filename);

  std::ostream &operator()(const level_t lvl);

 private:
  std::unique_ptr<stream> output_stream_;
  static std::size_t valid_rank_;
  static std::size_t rank_;
};

}
}
}

static parkway::utility::status::status_handler STATUS;
const auto info = parkway::utility::status::level_t::info;
const auto progress = parkway::utility::status::level_t::progress;

#endif  // UTILITY_LOGGING_HPP_
