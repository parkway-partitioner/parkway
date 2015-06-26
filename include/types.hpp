#ifndef TYPES_HPP_
#define TYPES_HPP_

namespace parkway {

// # visit order
enum class visit_order_t : uint_fast8_t {
  INCREASING = 1, DECREASING, RANDOM, INCREASING_WEIGHT, DECREASING_WEIGHT};

inline const char * to_c_str(const visit_order_t &item) {
  const char *str_array[] {"increasing", "decreasing", "random",
                           "increasing weight", "decreasing weight"};
  return str_array[static_cast<int>(item) - 1];
}


}

#endif  // TYPES_HPP_
