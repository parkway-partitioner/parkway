#ifndef UTILITY_PARKWAY_HPP_
#define UTILITY_PARKWAY_HPP_

namespace parkway {
namespace utility {

template <typename Type>
inline void quick_sort(const int left, const int right, Type *array) {
  int left_arrow = left;
  int right_arrow = right;
  Type pivot = array[(left + right) / 2];

  do {
    while (array[right_arrow] > pivot) {
      --right_arrow;
    }
    while (array[left_arrow] < pivot) {
      ++left_arrow;
    }

    if (left_arrow <= right_arrow) {
      std::swap(array[left_arrow], array[right_arrow]);
      ++left_arrow;
      --right_arrow;
    }
  } while (right_arrow >= left_arrow);

  if (left < right_arrow) {
    quick_sort(left, right_arrow, array);
  }

  if (left_arrow < right) {
    quick_sort(left_arrow, right, array);
  }
}

}
}


#endif  // UTILITY_PARKWAY_HPP_
