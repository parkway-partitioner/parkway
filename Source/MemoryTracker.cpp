#include <iostream>
#include <map>
#include "MemoryTracker.hpp"

using namespace std;

map<int, int> MemoryTracker::table;
int MemoryTracker::allocated = 0;
int MemoryTracker::max = 0;
int MemoryTracker::power = 0;

void *operator new(size_t n) throw(bad_alloc) {
	return MemoryTracker::allocate(n);
}

void *operator new[](size_t n) throw(bad_alloc) {
	return MemoryTracker::allocate(n);
}

void operator delete[](void *p) throw() {
	MemoryTracker::deallocate(p);
}

void operator delete(void *p) throw() {
	MemoryTracker::deallocate(p);
}

