#include "gtest/gtest.h"
#include "data_structures/map_to_pos_int.hpp"
#include "data_structures/internal/table_utils.hpp"

using parkway::data_structures::map_to_pos_int;

TEST(MapToPosInt, DefaultConstruction) {
  map_to_pos_int int_map;
  ASSERT_EQ(int_map.size(), 0);
}


TEST(MapToPosInt, Construction) {
  bool use_hash = false;
  int capacity = 10;
  map_to_pos_int int_map(capacity, use_hash);
  ASSERT_EQ(int_map.size(), 0);
  ASSERT_EQ(int_map.capacity(), capacity);
  ASSERT_EQ(int_map.use_hash(), use_hash);
}


TEST(MapToPosInt, InsertNoHash) {
  bool use_hash = false;
  int capacity = 3;
  map_to_pos_int int_map(capacity, use_hash);

  // False as no collision.
  ASSERT_FALSE(int_map.insert(0, 5));
  ASSERT_EQ(int_map.get(0), 5);
  // True as overwriting occurs.
  ASSERT_TRUE(int_map.insert(0, 4));
  ASSERT_EQ(int_map.get(0), 4);
}


TEST(MapToPosInt, InsertHash) {
  bool use_hash = true;
  int capacity = 8;
  map_to_pos_int int_map(capacity, use_hash);

  // False as no collision.
  ASSERT_FALSE(int_map.insert(0, 5));
  ASSERT_EQ(int_map.get(0), 5);
  // True as overwites old value (as no scatter array is set).
  ASSERT_TRUE(int_map.insert(0, 4));
  ASSERT_EQ(int_map.get(0), 4);

  parkway::data_structures::internal::table_utils::set_scatter_array(1);
  // Scatter is size 1 so the only item is index 0 -- overwrites index 0.
  ASSERT_TRUE(int_map.insert(0, 7));
  ASSERT_EQ(int_map.get(0), 7);
}


TEST(MapToPosInt, InsertIfEmptyNoHash) {
  bool use_hash = false;
  int capacity = 3;
  map_to_pos_int int_map(capacity, use_hash);

  // -1 if nothing in the slot.
  ASSERT_EQ(int_map.insert_if_empty(0, 5), -1);
  // Returns the item already in that slot.
  ASSERT_EQ(int_map.insert_if_empty(0, 4), 5);
}


TEST(MapToPosInt, Clear) {
  bool use_hash = false;
  int capacity = 3;
  map_to_pos_int int_map(capacity, use_hash);

  ASSERT_FALSE(int_map.insert(0, 5));
  ASSERT_FALSE(int_map.insert(1, 5));
  ASSERT_FALSE(int_map.insert(2, 5));
  ASSERT_EQ(int_map.size(), 3);
  int_map.clear();
  ASSERT_EQ(int_map.size(), 0);
}
