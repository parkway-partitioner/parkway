#ifndef UTILITY_STATUS_HPP_
#define UTILITY_STATUS_HPP_
#include <cstdarg>
#include <iostream>
#include <fstream>
#include <memory>

namespace parkway {
namespace utility {
namespace status {

class handler {
 public:
  static void set_rank(const std::size_t rank) {
    rank_ = rank;
  }

  static void set_rank_to_write_on(const std::size_t rank) {
    write_on_ = rank;
  }

  static void enable_info() {
    show_info_ = true;
  }

  static void disable_info() {
    show_info_ = false;
  }

  static bool info_enabled() {
    return show_info_;
  }

  static void enable_progress() {
    show_progress_ = true;
  }

  static void disable_progress() {
    show_progress_ = false;
  }

  static bool progress_enabled() {
    return show_progress_;
  }

  ~handler() {
    out_->flush();
    if (filestream_ != nullptr) {
      filestream_->close();
      delete filestream_;
    }
  }

  static void print_to_terminal() {
    if (filestream_ != nullptr) {
      filestream_->flush();
      filestream_->close();
      delete filestream_;
      filestream_ = nullptr;
    }
    out_ = &std::cout;
  }

  static bool set_log_file(const char *filename) {
    if (filestream_ != nullptr) {
      filestream_->close();
    }
    filestream_ = new std::ofstream(filename);
    if (filestream_->is_open()) {
      out_ = filestream_;
      return true;
    } else {
      return false;
    }
  }

  friend inline void info(const char *format, ...);
  friend inline void progress(const char *format, ...);
  friend inline void warning(const char *format, ...);
  friend inline void warning_on_processor(const char *format, ...);
  friend inline void error_on_processor(const char *format, ...);

 private:
  static std::ostream *out_;
  static std::ofstream *filestream_;
  static std::unique_ptr<char> buffer_;
  static std::size_t write_on_;
  static std::size_t rank_;
  static bool show_info_;
  static bool show_progress_;
};


inline void info(const char *format, ...) {
  if (handler::rank_ == handler::write_on_ &&
      handler::show_info_ &&
      handler::out_ != nullptr) {
    va_list args;
    va_start(args, format);
    std::vsprintf(handler::buffer_.get(), format, args);
    va_end(args);
    *handler::out_ << handler::buffer_.get();
    handler::out_->flush();
  }
}

inline void progress(const char *format, ...) {
  if (handler::rank_ == handler::write_on_ &&
      handler::show_progress_ &&
      handler::out_ != nullptr) {
    va_list args;
    va_start(args, format);
    std::vsprintf(handler::buffer_.get(), format, args);
    va_end(args);
    *handler::out_ << handler::buffer_.get();
    handler::out_->flush();
  }
}

inline void warning(const char *format, ...) {
  if (handler::rank_ == handler::write_on_ &&
      handler::show_info_ &&
      handler::out_ != nullptr) {
    char *buffer = handler::buffer_.get();
    std::sprintf(buffer, "[Warning!] ");
    va_list args;
    va_start(args, format);
    std::vsprintf(buffer, format, args);
    va_end(args);
    *handler::out_ << handler::buffer_.get();
  }
}

inline void warning_on_processor(const char *format, ...) {
  if (handler::out_ != nullptr) {
    char *buffer = handler::buffer_.get();
    std::sprintf(buffer, "[Warning!] ");
    va_list args;
    va_start(args, format);
    std::vsprintf(buffer, format, args);
    va_end(args);
    *handler::out_ << handler::buffer_.get();
  }
}

inline void error_on_processor(const char *format, ...) {
  if (handler::out_ != nullptr) {
    va_list args;
    char *buffer = handler::buffer_.get();
    std::sprintf(buffer, "\33[31m");
    va_start(args, format);
    std::vsprintf(buffer, format, args);
    va_end(args);
    std::sprintf(buffer, "\33[0m");
    *handler::out_ << handler::buffer_.get();
  }
}

}
}
}

#endif  // UTILITY_STATUS_HPP_
