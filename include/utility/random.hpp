#ifndef UTILITY_RANDOM_HPP_
#define UTILITY_RANDOM_HPP_

namespace parkway {
namespace utility {

inline int random(const int lower, const int upper) {
  if (lower == upper) {
    return lower;
  }
#ifdef USE_SPRNG
  return static_cast<int>(sprng() * (upper - lower) + lower);
#else
  return static_cast<int>(drand48() * (upper - lower) + lower);
#endif
}


template <typename Type>
inline void random_permutation(Type *array, int size) {
  for (int i = 0; i < size; ++i) {
    int j = random(i, size);
    std::swap(array[i], array[j]);
  }
}

}
}

#endif  // UTILITY_RANDOM_HPP_
