#include "gtest/gtest.h"
#include "utility/math.hpp"

namespace util = parkway::utility;

TEST(Math, Log2) {
  ASSERT_EQ(util::math::log2(0), 0);
  ASSERT_EQ(util::math::log2(1), 0);
  ASSERT_EQ(util::math::log2(2), 1);
  ASSERT_EQ(util::math::log2(3), 1);
  ASSERT_EQ(util::math::log2(16), 4);
}


TEST(Math, IsPowerOfTwo) {
  ASSERT_TRUE(util::math::is_power_of_2(0));
  ASSERT_TRUE(util::math::is_power_of_2(1));
  ASSERT_TRUE(util::math::is_power_of_2(2));
  ASSERT_TRUE(util::math::is_power_of_2(4));
  ASSERT_TRUE(util::math::is_power_of_2(1024));
  ASSERT_TRUE(util::math::is_power_of_2(65536));

  ASSERT_FALSE(util::math::is_power_of_2(-1));
  ASSERT_FALSE(util::math::is_power_of_2(3));
  ASSERT_FALSE(util::math::is_power_of_2(9));
  ASSERT_FALSE(util::math::is_power_of_2(1023));
}


TEST(Math, Half) {
  // Uses bit-shifting on ints.
  ASSERT_EQ(util::math::half(2), 1);
  ASSERT_EQ(util::math::half(3), 1);
  ASSERT_EQ(util::math::half(4), 2);
  ASSERT_EQ(util::math::half(1024), 512);
}
