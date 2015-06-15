#include <stdlib.h>
#include <iostream>
#include <cassert>
#include <map>

using namespace std;

class memory_tracker {
private:
  static map<int, int> table_;
  static int allocated_;
  static int maximum_;
  static int power_;

public:
  inline static void start() {
    power_ = 1;
    allocated_ = 0;
    maximum_ = 0;
  }

  inline static void stop() { power_ = 0; }

  inline static void add(void *p, int size) {
    int _power = power_;
    power_ = 0;
    table_[intptr_t(p)] = size;
    power_ = _power;
  }

  inline static int size(void *p) { return table_[intptr_t(p)]; }

  inline static void *allocate(size_t n) {
    void *p = malloc(n);
    assert(p);
    if (power_) {
      add(p, n);
      allocated_ += n;
      if (allocated_ > maximum_)
        maximum_ = allocated_;
    }
    return p;
  }

  inline static void deallocate(void *p) {
    if (!p)
      return;
    if (power_)
      allocated_ -= size(p);
    free(p);
  }

  inline static double usage() { return (double) allocated_ / 1024.0; }

  inline static double peak() { return (double) maximum_ / 1024.0; }
};
