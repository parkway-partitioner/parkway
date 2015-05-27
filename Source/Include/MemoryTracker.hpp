#include <stdlib.h>
#include <iostream>
#include <cassert>
#include <map>

using namespace std;

class MemoryTracker
{
private:

  static map<int,int> table;
  static int allocated;
  static int max;
  static int power;

public:

  inline static void start() {
    power = 1;
    allocated = 0;
    max = 0;
  }
  inline static void stop() {
    power = 0;
  }
  inline static void add(void *p, int size) {
    int _power = power;
    power = 0;
    table[*static_cast<int*>(p)] = size;
    power = _power;
  }
  inline static int size(void *p) {
    return table[*static_cast<int*>(p)];
  }
  inline static void *allocate(size_t n) {
    void *p = malloc(n);
    assert(p);
    if (power) {
      add(p,n);
      allocated += n;
      if (allocated > max)
	max = allocated;
    }
    return p;
  }
  inline static void deallocate(void *p) {
    if (!p)
      return;
    if (power)
      allocated -= size(p);
    free(p);
  }
  inline static double usage() {
    return (double) allocated/1024.0;
  }
  inline static double peak() {
    return (double) max/1024.0;
  }

};

void *operator new(size_t n);
void *operator new[](size_t n);
void operator delete[](void *p) throw();
void operator delete(void *p) throw();




