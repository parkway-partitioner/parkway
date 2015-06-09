#ifndef _DATA_STRUCTURES_MAP_FROM_POS_INT_HPP
#define _DATA_STRUCTURES_MAP_FROM_POS_INT_HPP

#include "data_structures/dynamic_array.hpp"
#include "data_structures/internal/table_utils.hpp"

namespace parkway {
namespace data_structures {

/* New MapFromPosInt Class */
/* requires positive int keys */
template<typename T> class map_from_pos_int {
 protected:
  int numEntries;
  int size;

  dynamic_array<T> table;
  dynamic_array<int> keys;

 public:
  map_from_pos_int() : map_from_pos_int(0) {
  }

  map_from_pos_int(int _size) {
    createTable(_size);
  }

  ~map_from_pos_int() {
  }

  void createTable(int _size) {
    numEntries = 0;
    size = internal::table_utils::table_size(_size);

    #ifdef DEBUG_TABLES
    assert(size >= _size);
    #endif

    keys.reserve(size);
    table.reserve(size);

    for (std::size_t i = 0; i < size; ++i) {
      keys[i] = -1;
    }
  }

  void destroyTable() {
    numEntries = 0;
    table.reserve(0);
    keys.reserve(0);
  }

  void recoverTable() {
    table.reserve(size);
    keys.reserve(size);
    numEntries = 0;

    for (std::size_t i = 0; i < size; ++i) {
      keys[i] = -1;
    }
  }

  int insertKey(int key, T val) {
    int indepKey = internal::table_utils::scatter_key(key);
    int slot = internal::hashes::primary(indepKey, size);

    while (keys[slot] != -1 && keys[slot] != indepKey)
      slot = (slot + internal::hashes::secondary(indepKey, size)) % size;

    if (keys[slot] == -1) {
      table[slot] = val;
      keys[slot] = indepKey;
      ++numEntries;
      return 0;
    } else {
      table[slot] = val;
      return 1;
    }
  }


  T &getVal(int key) {
    int indepKey = internal::table_utils::scatter_key(key);
    int slot = internal::hashes::primary(indepKey, size);
    while (keys[slot] != indepKey) {
      slot = (slot + internal::hashes::secondary(indepKey, size)) % size;
      #ifdef DEBUG_TABLES
      assert(keys[slot] != -1);
      #endif
    }

    return table[slot];
  }

  inline int getNumEntries() {
    return numEntries;
  }

  inline int getNumSlots() {
    return size;
  }
};

}  // namespace data_structures
}  // namespace parkway

#endif  // _DATA_STRUCTURES_MAP_FROM_POS_INT_HPP
