#include "gtest/gtest.h"
#include "data_structures/bit_field.hpp"

using parkway::data_structures::bit_field;

TEST(BitField, Construction) {
  bit_field bf1;

  ASSERT_EQ(bf1.capacity(), 0);
  ASSERT_EQ(bf1.number_of_bits(), 0);

  bit_field bf2(10);
  ASSERT_EQ(bf2.capacity(), 1);
  ASSERT_EQ(bf2.number_of_bits(), 10);
}


TEST(BitField, OperatorBrackets) {
  bit_field bf(1);
  ASSERT_FALSE(bf(0));
  ASSERT_FALSE(bf[0]);
}


TEST(BitField, Set) {
  bit_field bf(4);
  ASSERT_FALSE(bf[2]);
  bf.set(2);
  ASSERT_TRUE(bf[2]);
}


TEST(BitField, Unset) {
  bit_field bf(4);
  bf.set(3);
  ASSERT_TRUE(bf[3]);
  bf.unset(3);
  ASSERT_FALSE(bf[3]);
}


TEST(BitField, TestAll) {
  bit_field bf(10);
  ASSERT_FALSE(bf.test_all());

  bf.set();
  ASSERT_TRUE(bf.test_all());
}


TEST(BitField, Check) {
  bit_field bf;
  ASSERT_EQ(bf.capacity(), 0);
  // Check if there is space in index 0.
  bf.check(0);
  // Should now space for bits.
  ASSERT_NE(bf.capacity(), 0);
}


TEST(BitField, Chunk) {
  bit_field bf(64);
  // Set all bits.
  bf.set();
  ASSERT_EQ(bf.chunk(0), bf.CHUNK_MAX);
  ASSERT_EQ(bf.chunk(1), bf.CHUNK_MAX);
}


TEST(BitField, Data) {
  bit_field bf(64);
  bf.set();

  bit_field::chunk_t *chunks = bf.data();
  chunks[0] = 0;

  ASSERT_EQ(bf.chunk(0), 0);
  ASSERT_EQ(bf.chunk(1), bf.CHUNK_MAX);
}


TEST(BitField, Reserve) {
  bit_field bf;
  ASSERT_EQ(bf.capacity(), 0);
  ASSERT_EQ(bf.number_of_bits(), 0);

  // Reserve 1 less than the width of the chunk, so only 1 chunk will be
  // reserved (i.e. capacity == 1).
  bf.reserve(bf.CHUNK_WIDTH - 1);
  ASSERT_EQ(bf.capacity(), 1);
  ASSERT_EQ(bf.number_of_bits(), bf.CHUNK_WIDTH - 1);
}


TEST(BitField, SetUnsetAll) {
  bit_field bf(64);
  ASSERT_FALSE(bf.test_all());

  bf.set();
  ASSERT_TRUE(bf.test_all());

  bf.unset();
  ASSERT_FALSE(bf.test_all());
}
