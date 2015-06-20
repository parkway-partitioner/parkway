#include "gtest/gtest.h"
#include "data_structures/dynamic_array.hpp"

using parkway::data_structures::dynamic_array;

TEST(DynamicArray, Construction) {
  dynamic_array<int> arr1;
  ASSERT_EQ(arr1.capacity(), 0);

  dynamic_array<int> arr2(5);
  ASSERT_EQ(arr2.capacity(), 5);
}


TEST(DynamicArray, Reserve) {
  dynamic_array<double> arr;
  ASSERT_EQ(arr.capacity(), 0);
  arr.reserve(0);
  ASSERT_EQ(arr.capacity(), 0);
  arr.reserve(5);
  ASSERT_EQ(arr.capacity(), 5);
}



TEST(DynamicArray, Capacity) {
  dynamic_array<int> arr(12);
  ASSERT_EQ(arr.capacity(), 12);
}


TEST(DynamicArray, OperatorBracketAccess) {
  dynamic_array<int> arr(5);
  arr[0] = 10;
  ASSERT_EQ(arr[0], 10);
}


// TODO(gb610): think about memory usage, this causes the tests to break as
// the pointer to arr1's underlying data is copied to arr2. When arr1 and arr2
// are deleted this memory is freed, when the other is deleted an already freed
// array array tries to get deleteed.
//
// TEST(DynamicArray, SetData) {
//   dynamic_array<int> arr1(5);
//   arr1.assign(0, 10);
//   arr1.assign(1, 7);
//   arr1.assign(2, 6);
//
//   dynamic_array<int> arr2;
//   arr2.set_data(arr1.data(), arr1.capacity());
//
//   ASSERT_EQ(arr2.data(), arr1.data());
//   ASSERT_EQ(arr2.capacity(), arr1.capacity());
// }

