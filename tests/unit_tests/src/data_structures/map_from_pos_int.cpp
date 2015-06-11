#include "gtest/gtest.h"
#include "data_structures/map_from_pos_int.hpp"
#include "data_structures/internal/table_utils.hpp"

using parkway::data_structures::map_from_pos_int;
using parkway::data_structures::internal::table_utils;

TEST(MapFromPosInt, DefaultConstruction) {
  map_from_pos_int<double> map_;
  ASSERT_EQ(map_.size(), 0);
  ASSERT_GE(map_.capacity(), 0);
}

TEST(MapFromPosInt, ConstructWithSize) {
  map_from_pos_int<char> map_(10);
  ASSERT_EQ(map_.size(), 0);
  ASSERT_GE(map_.capacity(), 10);
}

TEST(MapFromPosInt, Create) {
  map_from_pos_int<float> map_;
  int old_capacity = map_.capacity();
  int new_capacity = old_capacity + 1;
  map_.create(new_capacity);
  ASSERT_GT(map_.capacity(), old_capacity);
  ASSERT_GE(map_.capacity(), new_capacity);
}

TEST(MapFromPosInt, Insert) {
  map_from_pos_int<char> map_(10);
  table_utils::set_scatter_array(10);

  ASSERT_EQ(map_.size(), 0);
  ASSERT_FALSE(map_.insert(0, 'a'));
  ASSERT_EQ(map_.size(), 1);
  // True as a value is overwritten.
  ASSERT_TRUE(map_.insert(0, 'b'));
  ASSERT_EQ(map_.size(), 1);

  ASSERT_FALSE(map_.insert(5, 'c'));
  ASSERT_EQ(map_.size(), 2);

  // Unset scatter key, other tests rely on it being unset.
  table_utils::set_scatter_array(table_utils::SCATTER_ARRAY_NOT_SET);
}

TEST(MapFromPosInt, Get) {
  map_from_pos_int<char> map_(10);
  table_utils::set_scatter_array(10);

  ASSERT_EQ(map_.size(), 0);
  ASSERT_FALSE(map_.insert(0, 'a'));
  ASSERT_EQ(map_.size(), 1);
  ASSERT_FALSE(map_.insert(5, 'c'));
  ASSERT_EQ(map_.size(), 2);

  ASSERT_EQ(map_.get(0), 'a');
  ASSERT_EQ(map_[5], 'c');

  // Unset scatter key, other tests rely on it being unset.
  table_utils::set_scatter_array(table_utils::SCATTER_ARRAY_NOT_SET);
}


TEST(MapFromPosInt, Destroy) {
  map_from_pos_int<float> map_;
  table_utils::set_scatter_array(map_.capacity());

  ASSERT_FALSE(map_.insert(0, 1.23));
  ASSERT_FALSE(map_.insert(1, 4.56));
  ASSERT_FALSE(map_.insert(2, 7.89));
  ASSERT_EQ(map_.size(), 3);

  map_.destroy();
  ASSERT_EQ(map_.size(), 0);

  // Unset scatter key, other tests rely on it being unset.
  table_utils::set_scatter_array(table_utils::SCATTER_ARRAY_NOT_SET);
}
