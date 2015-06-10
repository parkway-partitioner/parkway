#include "gtest/gtest.h"
#include "data_structures/dynamic_array.hpp"

using parkway::data_structures::dynamic_array;

TEST(DynamicArray, Construction) {
  dynamic_array<int> arr1;
  ASSERT_EQ(arr1.capacity(), 0);
  ASSERT_EQ(arr1.data(), nullptr);

  dynamic_array<int> arr2(-5);
  ASSERT_EQ(arr2.capacity(), 0);
  ASSERT_EQ(arr2.data(), nullptr);

  dynamic_array<int> arr3(5);
  ASSERT_EQ(arr3.capacity(), 5);
  ASSERT_NE(arr3.data(), nullptr);
}


TEST(DynamicArray, Reserve) {
  dynamic_array<double> arr;
  ASSERT_EQ(arr.capacity(), 0);
  ASSERT_TRUE(arr.reserve(0));
  ASSERT_EQ(arr.capacity(), 0);

  ASSERT_TRUE(arr.reserve(5));
  ASSERT_EQ(arr.capacity(), 5);
}


TEST(DynamicArray, Check) {
  dynamic_array<int> arr;
  int capacity = arr.capacity();
  // Checking an out of range index should expand the array.
  arr.check(capacity + 1);
  ASSERT_NE(arr.capacity(), capacity);
}


TEST(DynamicArray, Capacity) {
  dynamic_array<int> arr(12);
  ASSERT_EQ(arr.capacity(), 12);
}


TEST(DynamicArray, Assign) {
  dynamic_array<int> arr(5);
  arr.assign(0, 10);
  ASSERT_EQ(arr(0), 10);
}

TEST(DynamicArray, OperatorBracketAccess) {
  dynamic_array<int> arr(5);
  arr.assign(0, 10);
  ASSERT_EQ(arr(0), 10);
  ASSERT_EQ(arr[0], 10);
}


TEST(DynamicArray, Find) {
  dynamic_array<int> arr(5);
  arr.assign(0, 10);
  arr.assign(1, 7);
  ASSERT_EQ(arr.find(7), 1);
  ASSERT_EQ(arr.find(10), 0);

  ASSERT_EQ(arr.find(5), arr.not_found());
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


TEST(DynamicArray, Adjust) {
  dynamic_array<int> arr;
  ASSERT_EQ(arr.capacity(), 0);
  ASSERT_EQ(arr.data(), nullptr);

  arr.adjust(-10);
  ASSERT_EQ(arr.capacity(), 0);
  ASSERT_EQ(arr.data(), nullptr);

  arr.adjust(0);
  ASSERT_EQ(arr.capacity(), 0);
  ASSERT_EQ(arr.data(), nullptr);

  arr.adjust(10);
  ASSERT_EQ(arr.capacity(), 10);
  ASSERT_NE(arr.data(), nullptr);

  arr.adjust(-10);
  ASSERT_EQ(arr.capacity(), 0);
  ASSERT_EQ(arr.data(), nullptr);
}



