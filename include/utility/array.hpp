#ifndef UTILITY_ARRAY_HPP_
#define UTILITY_ARRAY_HPP_

#include "data_structures/dynamic_array.hpp"

namespace parkway {
namespace utility {
namespace ds = parkway::data_structures;

template <typename T>
inline void set_to(T* data, const std::size_t length, const T value) {
  std::fill(data, data + length, value);
}

template <typename T>
inline void set_to_zero(T* data, const std::size_t length) {
  set_to(data, length, 0);
}

}  // namespace utility
}  // namespace parkway

#endif  // UTILITY_ARRAY_HPP_
