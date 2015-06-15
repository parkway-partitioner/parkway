#include <iostream>
#include <map>
#include <new>
#include "memory_tracker.hpp"

using namespace std;

map<int, int> memory_tracker::table_;
int memory_tracker::allocated_ = 0;
int memory_tracker::maximum_ = 0;
int memory_tracker::power_ = 0;

void *operator new(std::size_t n) {
  return memory_tracker::allocate(n);
}

void *operator new[](std::size_t n) {
  return memory_tracker::allocate(n);
}

void operator delete[](void *p) throw() { memory_tracker::deallocate(p); }

void operator delete(void *p) throw() { memory_tracker::deallocate(p); }
