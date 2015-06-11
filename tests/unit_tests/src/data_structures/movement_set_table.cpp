#include "gtest/gtest.h"
#include "data_structures/movement_set_table.hpp"

using parkway::data_structures::movement_set_table;
using parkway::data_structures::movement_set;

TEST(MovementSet, Construction) {
  movement_set ms(1, 2, 3);
  ASSERT_EQ(ms.gain, 1);
  ASSERT_EQ(ms.weight, 2);
  ASSERT_EQ(ms.proc, 3);
}

TEST(MovementSetTable, Construction) {
  movement_set_table mst(8, 4);
  ASSERT_EQ(mst.number_of_parts(), 8);
  ASSERT_EQ(mst.number_of_processors(), 4);
}

TEST(MovementSetTable, SetMaxPartWeight) {
  movement_set_table mst(1, 2);
  mst.set_max_part_weight(20);
  ASSERT_EQ(mst.max_part_weight(), 20);
}
