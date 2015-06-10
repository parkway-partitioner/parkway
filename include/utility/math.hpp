#ifndef _UTILITY_MATH_HPP
#define _UTILITY_MATH_HPP

namespace parkway {
namespace utility {
namespace math {

inline static int log2(int of) {
  int result = 0;
  while (of >>= 1) {
    ++result;
  }
  return result;
}


inline static bool is_power_of_2(int num) {
  if (num == 1) {
    return true;
  } else if (num < 0) {
    return false;
  }

  while (num > 1) {
    // If any zeros in binary representation then cannot be a power of 2.
    if (num & 0x1) {
      return false;
    }
    num >>= 1;
  }

  return true;
}


inline int half(int number) {
  return number >> 1;
}

}  // namespace math
}  // namespace utility
}  // namespace parkway

#endif  // _UTILITY_MATH_HPP
