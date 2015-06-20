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

enum class sort_order {INCREASING, DECREASING};

template <typename Type>
inline void quick_sort_by_another_array(const int left, const int right,
                                        Type *array, const Type *val_array,
                                        sort_order order) {
  int left_arrow = left;
  int right_arrow = right;
  int pivot = array[(left + right) / 2];

  if (order == sort_order::INCREASING) {
    do {
      while (val_array[array[right_arrow]] > val_array[pivot]) {
        --right_arrow;
      }

      while (val_array[array[left_arrow]] < val_array[pivot]) {
        ++left_arrow;
      }

      if (left_arrow <= right_arrow) {
        std::swap(array[left_arrow++], array[right_arrow--]);
      }

    } while (right_arrow >= left_arrow);

    if (left < right_arrow) {
      quick_sort_by_another_array(left, right_arrow, array, val_array, order);
    }

    if (left_arrow < right) {
      quick_sort_by_another_array(left_arrow, right, array, val_array, order);
    }

  } else {
    do {
      while (val_array[array[right_arrow]] < val_array[pivot]) {
        --right_arrow;
      }
      while (val_array[array[left_arrow]] > val_array[pivot]) {
        ++left_arrow;
      }

      if (left_arrow <= right_arrow) {
        std::swap(array[left_arrow++], array[right_arrow--]);
      }
    } while (right_arrow >= left_arrow);

    if (left < right_arrow) {
      quick_sort_by_another_array(left, right_arrow, array, val_array, order);
    }

    if (left_arrow < right) {
      quick_sort_by_another_array(left_arrow, right, array, val_array, order);
    }
  }
}


}
}


#endif  // UTILITY_PARKWAY_HPP_
