#include "gtest/gtest.h"
#include "data_structures/bit_field.hpp"

using parkway::data_structures::bit_field;

TEST(bit_field, construction) {
  bit_field bf1;

  ASSERT_EQ(bf1.capacity(), 0);
  ASSERT_EQ(bf1.number_of_bits(), 0);

  bit_field bf2(10);
  ASSERT_EQ(bf2.capacity(), 1);
  ASSERT_EQ(bf2.number_of_bits(), 10);
}


TEST(bit_field, subscript_operators) {
  bit_field bf(1);
  ASSERT_FALSE(bf(0));
  ASSERT_FALSE(bf[0]);
}


TEST(bit_field, set) {
  bit_field bf(4);
  ASSERT_FALSE(bf[2]);
  bf.set(2);
  ASSERT_TRUE(bf[2]);
}


TEST(bit_field, unset) {
  bit_field bf(4);
  bf.set(3);
  ASSERT_TRUE(bf[3]);
  bf.unset(3);
  ASSERT_FALSE(bf[3]);
}


TEST(bit_field, test_all) {
  bit_field bf(10);
  ASSERT_FALSE(bf.test_all());

  bf.set();
  ASSERT_TRUE(bf.test_all());
}


TEST(bit_field, check) {
  bit_field bf;
  ASSERT_EQ(bf.capacity(), 0);
  // Check if there is space in index 0.
  bf.check(0);
  // Should now space for bits.
  ASSERT_NE(bf.capacity(), 0);
}


TEST(bit_field, chunk) {
  bit_field bf(64);
  // Set all bits.
  bf.set();
  ASSERT_EQ(bf.chunk(0), bf.CHUNK_MAX);
  ASSERT_EQ(bf.chunk(1), bf.CHUNK_MAX);
}


TEST(bit_field, data) {
  bit_field bf(64);
  bf.set();

  auto chunks = bf.data();
  chunks[0] = 0;

  ASSERT_EQ(bf.chunk(0), 0);
  ASSERT_EQ(bf.chunk(1), bf.CHUNK_MAX);
}


TEST(bit_field, reserve) {
  bit_field bf;
  ASSERT_EQ(bf.capacity(), 0);
  ASSERT_EQ(bf.number_of_bits(), 0);

  // Reserve 1 less than the width of the chunk, so only 1 chunk will be
  // reserved (i.e. capacity == 1).
  bf.reserve(bf.CHUNK_WIDTH - 1);
  ASSERT_EQ(bf.capacity(), 1);
  ASSERT_EQ(bf.number_of_bits(), bf.CHUNK_WIDTH - 1);
}


TEST(bit_field, set_all) {
  bit_field bf(64);
  ASSERT_FALSE(bf.test_all());

  bf.set();
  ASSERT_TRUE(bf.test_all());
}

TEST(bit_field, unset_all) {
  bit_field bf(64);
  bf.set();
  ASSERT_TRUE(bf.test_all());

  bf.unset();
  ASSERT_FALSE(bf.test_all());
}
