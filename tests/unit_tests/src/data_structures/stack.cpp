#include "gtest/gtest.h"
#include "data_structures/stack.hpp"

using parkway::data_structures::stack;

TEST(Stack, Size) {
  stack<int> int_stack;
  ASSERT_EQ(int_stack.size(), 0);
}


TEST(Stack, Push) {
  stack<int> int_stack;
  int_stack.push(4);
  ASSERT_EQ(int_stack.size(), 1);
}


TEST(Stack, Pop) {
  stack<int> int_stack;
  int_stack.push(4);
  ASSERT_EQ(int_stack.size(), 1);
  ASSERT_EQ(int_stack.pop(), 4);
  ASSERT_EQ(int_stack.size(), 0);
}


TEST(Stack, Top) {
  stack<int> int_stack;
  int_stack.push(4);
  ASSERT_EQ(int_stack.size(), 1);
  ASSERT_EQ(int_stack.top(), 4);
  ASSERT_EQ(int_stack.size(), 1);
}
