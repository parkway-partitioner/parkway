#include "gtest/gtest.h"
#include "data_structures/match_request_table.hpp"

using parkway::data_structures::match_request_table;

TEST(MatchRequestTable, Construction) {
  match_request_table table_(10);
  ASSERT_EQ(table_.size(), 0);
  ASSERT_EQ(table_.capacity(), 10);
}
