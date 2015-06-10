#include "gtest/gtest.h"
#include "data_structures/complete_binary_tree.hpp"

using parkway::data_structures::complete_binary_tree;

TEST(CompleteBinaryTree, EmptyConstruction) {
  complete_binary_tree<int> bt;
  ASSERT_EQ(bt.size(), 0);
}


TEST(CompleteBinaryTree, Construction) {
  int values[] = {3, 4, 5};
  int keys[] = {0, 1, 2};
  int entries = 3;

  complete_binary_tree<int> bt(values, keys, entries);
  ASSERT_EQ(bt.size(), entries);
}


TEST(CompleteBinaryTree, Setup) {
  complete_binary_tree<int> bt;
  ASSERT_EQ(bt.size(), 0);

  int values[] = {6, 7, 8, 9};
  int keys[] = {1, 2, 3, 4};
  int entries = 4;

  bt.setup(values, keys, entries);
  ASSERT_EQ(bt.size(), entries);
}
