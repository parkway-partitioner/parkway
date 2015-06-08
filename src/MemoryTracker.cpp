#include <iostream>
#include <map>
#include <new>
#include "MemoryTracker.hpp"

using namespace std;

map<int, int> MemoryTracker::table;
int MemoryTracker::allocated = 0;
int MemoryTracker::max = 0;
int MemoryTracker::power = 0;

void *operator new(std::size_t n) {
  return MemoryTracker::allocate(n);
}

void *operator new[](std::size_t n) {
  return MemoryTracker::allocate(n);
}

void operator delete[](void *p) throw() { MemoryTracker::deallocate(p); }

void operator delete(void *p) throw() { MemoryTracker::deallocate(p); }
