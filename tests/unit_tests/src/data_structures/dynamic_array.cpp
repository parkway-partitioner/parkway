#include "gtest/gtest.h"
#include "data_structures/dynamic_array.hpp"

using parkway::data_structures::dynamic_array;

TEST(dynamic_array, construction) {
  dynamic_array<int> arr1;
  ASSERT_EQ(arr1.capacity(), 0);

  dynamic_array<int> arr2(5);
  ASSERT_EQ(arr2.capacity(), 5);

  dynamic_array<int> arr3(arr1);
  ASSERT_EQ(arr1.size(), arr3.size());
  ASSERT_EQ(arr1.capacity(), arr3.capacity());
  ASSERT_EQ(arr1.data(), arr3.data());
}


TEST(dynamic_array, assignment_operator_by_reference) {
  dynamic_array<int> arr1;
  dynamic_array<int> arr2 = arr1;
  ASSERT_EQ(arr1.data(), arr2.data());
}

TEST(dynamic_array, assignment_operator_by_move) {
  dynamic_array<int> arr1 = dynamic_array<int>(3, 4);
  ASSERT_EQ(arr1.size(), 3);
  for (std::size_t i = 0; i < 3; ++i) {
    ASSERT_EQ(arr1[i], 4);
  }
}

TEST(dynamic_array, assign_sets_all_values) {
  int size = 5;
  int value = 10;
  dynamic_array<int> arr1;
  arr1.assign(size, value);
  ASSERT_EQ(arr1.size(), size);
  for (std::size_t i = 0; i < arr1.size(); ++i) {
    ASSERT_EQ(arr1[i], value);
  }
}

TEST(dynamic_array, assign_changes_size) {
  int size = 5;
  int init_value = 0;
  int value = 10;
  dynamic_array<int> arr1(2 * size, init_value);
  ASSERT_EQ(arr1.size(), 2 * size);
  for (std::size_t i = 0; i < arr1.size(); ++i) {
    ASSERT_EQ(arr1[i], init_value);
  }

  arr1.assign(size, value);
  ASSERT_EQ(arr1.size(), size);
  for (std::size_t i = 0; i < arr1.size(); ++i) {
    ASSERT_EQ(arr1[i], value);
  }
}


TEST(dynamic_array, at_resizes_array_if_too_small_in_assignment) {
  dynamic_array<int> arr;
  ASSERT_EQ(arr.size(), 0);
  ASSERT_EQ(arr.capacity(), 0);
  arr.at(50) = 10;
  ASSERT_EQ(arr.size(), 51);
  ASSERT_GE(arr.capacity(), 51);
  ASSERT_EQ(arr.at(50), 10);
}

TEST(dynamic_array, at_retrieves_value) {
  dynamic_array<int> arr(5, 3);
  ASSERT_EQ(arr.at(2), 3);
}

TEST(dynamic_array, subscript_resizes_array_if_too_small_in_assignment) {
  dynamic_array<int> arr;
  ASSERT_EQ(arr.size(), 0);
  ASSERT_EQ(arr.capacity(), 0);
  arr[50] = 10;
  ASSERT_EQ(arr.size(), 51);
  ASSERT_GE(arr.capacity(), 51);
  ASSERT_EQ(arr[50], 10);
}

TEST(dynamic_array, subscript_retrieves_value) {
  dynamic_array<int> arr(5, 3);
  ASSERT_EQ(arr[2], 3);
}

TEST(dynamic_array, front_sets_and_retrieves_value) {
  dynamic_array<int> arr(5, 3);
  arr.front() = 10;
  ASSERT_EQ(arr.front(), 10);
  ASSERT_EQ(&arr.front(), &arr[0]);
}

TEST(dynamic_array, back_sets_and_retrieves_value) {
  dynamic_array<int> arr(5, 3);
  arr.back() = 10;
  ASSERT_EQ(arr.back(), 10);
  ASSERT_EQ(&arr.back(), &arr[4]);
}

TEST(dynamic_array, data_retrieves_underlying_array) {
  dynamic_array<int> arr(5);
  for (std::size_t i = 0; i < arr.size(); ++i) {
    arr[i] = i;
  }
  int *data = arr.data();
  for (std::size_t i = 0; i < arr.size(); ++i) {
    ASSERT_EQ(&data[i], &arr[i]);
  }
}

TEST(dynamic_array, iterators_work_for_range_based_for_loop) {
  dynamic_array<int> arr(5);
  int i = 0;
  for (const int &item : arr) {
    ASSERT_EQ(&item, &arr[i++]);
  }
}

TEST(dynamic_array, empty_on_default_construction) {
  dynamic_array<int> arr;
  ASSERT_TRUE(arr.empty());
  arr.assign(5, 10);
  ASSERT_FALSE(arr.empty());
}

TEST(dynamic_array, reserve_can_increase_capacity) {
  dynamic_array<int> arr;
  ASSERT_EQ(arr.capacity(), 0);
  arr.reserve(50);
  ASSERT_EQ(arr.capacity(), 50);
}

TEST(dynamic_array, reserve_below_current_size_has_no_effect) {
  dynamic_array<int> arr(50);
  ASSERT_EQ(arr.size(), 50);
  ASSERT_EQ(arr.capacity(), 50);
  arr.reserve(0);
  ASSERT_EQ(arr.size(), 50);
  ASSERT_EQ(arr.capacity(), 50);
}

TEST(dynamic_array, max_size) {
  dynamic_array<int> arr;
  ASSERT_GT(arr.max_size(), 0);
}

TEST(dynamic_array, shrink_to_fit) {
  dynamic_array<int> arr;
  arr.reserve(50);
  ASSERT_EQ(arr.capacity(), 50);
  ASSERT_EQ(arr.size(), 0);
  arr.shrink_to_fit();
  ASSERT_EQ(arr.capacity(), 0);
  ASSERT_EQ(arr.size(), 0);
}

TEST(dynamic_array, clear_removes_all_elements) {
  dynamic_array<int> arr(10, 4);
  ASSERT_EQ(arr.capacity(), 10);
  ASSERT_EQ(arr.size(), 10);
  arr.clear();
  ASSERT_EQ(arr.capacity(), 10);
  ASSERT_EQ(arr.size(), 0);
}

TEST(dynamic_array, clear_and_shrink_reduces_capacity_to_zero) {
  dynamic_array<int> arr(10, 4);
  ASSERT_EQ(arr.capacity(), 10);
  ASSERT_EQ(arr.size(), 10);
  arr.clear_and_shrink();
  ASSERT_EQ(arr.capacity(), 0);
  ASSERT_EQ(arr.size(), 0);
}

TEST(dynamic_array, insert_const_reference_at_position) {
  dynamic_array<int> arr(10);
  ASSERT_EQ(arr.size(), 10);
  int new_value = 5;
  arr.insert(arr.end(), new_value);
  ASSERT_EQ(arr.size(), 11);
  ASSERT_EQ(arr[10], new_value);
}

TEST(dynamic_array, insert_move_at_position) {
  dynamic_array<int> arr(10);
  ASSERT_EQ(arr.size(), 10);
  arr.insert(arr.end(), 5);
  ASSERT_EQ(arr.size(), 11);
  ASSERT_EQ(arr[10], 5);
}

TEST(dynamic_array, set_data) {
  dynamic_array<int> arr;
  ASSERT_EQ(arr.size(), 0);
  ASSERT_EQ(arr.capacity(), 0);

  int new_data[] = {1, 2, 3, 4, 5};
  arr.set_data(new_data, 5);
  ASSERT_EQ(arr.size(), 5);
  for (std::size_t i = 0; i < arr.size(); ++i) {
    ASSERT_EQ(arr[i], new_data[i]);
  }
}

TEST(dynamic_array, erase_index_removes_an_element) {
  dynamic_array<int> arr(5);
  for (std::size_t i = 0; i < arr.size(); ++i) {
    arr[i] = i;
  }

  auto iter = arr.begin();
  iter += 2;
  arr.erase(iter);
  ASSERT_EQ(arr.size(), 4);
  ASSERT_EQ(arr[2], 3);
  ASSERT_EQ(arr[3], 4);
}

TEST(dynamic_array, erase_range_removes_multiple_elements) {
  dynamic_array<int> arr(5);
  for (std::size_t i = 0; i < arr.size(); ++i) {
    arr[i] = i;
  }

  auto iter = arr.begin();
  ++iter;
  arr.erase(iter, arr.end());
  ASSERT_EQ(arr.size(), 1);
  ASSERT_EQ(arr[0], 0);
}


TEST(dynamic_array, push_back_const_reference) {
  dynamic_array<int> arr;
  int value1 = 5;
  int value2 = 100;
  arr.push_back(value1);
  arr.push_back(value2);
  ASSERT_EQ(arr[0], 5);
  ASSERT_EQ(arr[1], 100);
  ASSERT_EQ(arr.size(), 2);
}


TEST(dynamic_array, push_back_move) {
  dynamic_array<int> arr;
  arr.push_back(5);
  arr.push_back(100);
  ASSERT_EQ(arr[0], 5);
  ASSERT_EQ(arr[1], 100);
  ASSERT_EQ(arr.size(), 2);
}

TEST(dynamic_array, pop_back_decreases_size) {
  dynamic_array<int> arr;
  arr.push_back(5);
  arr.push_back(1);
  ASSERT_EQ(arr.size(), 2);
  arr.pop_back();
  ASSERT_EQ(arr.size(), 1);
}

TEST(dynamic_array, sort_whole_array) {
  dynamic_array<int> arr;
  arr.push_back(1);
  arr.push_back(3);
  arr.push_back(4);
  arr.push_back(2);
  arr.push_back(0);
  ASSERT_EQ(arr.size(), 5);
  arr.sort();
  for (std::size_t i = 0; i < arr.size(); ++i) {
    EXPECT_EQ(arr[i], i);
  }
}

TEST(dynamic_array, sort_between_indices) {
  dynamic_array<int> arr;
  arr.push_back(1);
  arr.push_back(4);
  arr.push_back(3);
  arr.push_back(2);
  arr.push_back(0);
  ASSERT_EQ(arr.size(), 5);
  // First three elements -> {1, 3, 4, 2, 0}
  arr.sort_between(0, 2);
  ASSERT_EQ(arr[0], 1);
  ASSERT_EQ(arr[1], 3);
  ASSERT_EQ(arr[2], 4);

  // Sort while array.
  arr.sort_between(0, arr.size() - 1);
  for (std::size_t i = 0; i < arr.size(); ++i) {
    ASSERT_EQ(arr[i], i);
  }
}

TEST(dynamic_array, sort_using_another_array) {
  dynamic_array<int> arr1;
  for (std::size_t i = 0; i < 5; ++i) {
    arr1.push_back(i);
  }

  dynamic_array<int> sort_by;
  for (std::size_t i = 0; i < 5; ++i) {
    sort_by.push_back(4 - i);
  }

  arr1.sort_using_another_array(sort_by);
  for (std::size_t i = 0; i < 5; ++i) {
    ASSERT_EQ(arr1[i], 4 - i);
  }
}

TEST(dynamic_array, random_permutation) {
  dynamic_array<int> arr;
  for (std::size_t i = 0; i < 20; ++i) {
    arr.push_back(i);
  }
  arr.random_permutation();
  bool changed = false;
  for (std::size_t i = 0; i < arr.size(); ++i) {
    changed = changed || (arr[i] != i);
  }
  ASSERT_TRUE(changed);
}

TEST(dynamic_array, shared_ptr_works_on_assignment) {
  dynamic_array<int> arr;
  for (std::size_t i = 0; i < 20; ++i) {
    arr.push_back(i);
  }

  dynamic_array<int>* dynamic_arr = new dynamic_array<int>;
  *dynamic_arr = arr;

  ASSERT_EQ(dynamic_arr->size(), arr.size());
  ASSERT_EQ(dynamic_arr->capacity(), arr.capacity());
  ASSERT_EQ(dynamic_arr->data(), arr.data());
  int *underlying_data = dynamic_arr->data();

  delete dynamic_arr;

  ASSERT_EQ(arr.data(), underlying_data);
  ASSERT_NE(arr.data(), nullptr);
  ASSERT_EQ(arr.size(), 20);
}

TEST(dynamic_array, shared_ptr_works_on_copy_construction) {
  dynamic_array<int> arr;
  for (std::size_t i = 0; i < 20; ++i) {
    arr.push_back(i);
  }

  dynamic_array<int>* dynamic_arr = new dynamic_array<int>(arr);

  ASSERT_EQ(dynamic_arr->size(), arr.size());
  ASSERT_EQ(dynamic_arr->capacity(), arr.capacity());
  ASSERT_EQ(dynamic_arr->data(), arr.data());
  int *underlying_data = dynamic_arr->data();

  delete dynamic_arr;

  ASSERT_EQ(arr.data(), underlying_data);
  ASSERT_NE(arr.data(), nullptr);
  ASSERT_EQ(arr.size(), 20);
}
