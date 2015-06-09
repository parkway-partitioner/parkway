#ifndef _DATA_STRUCTURES_MAP_TO_POS_INT_HPP
#define _DATA_STRUCTURES_MAP_TO_POS_INT_HPP

#include "data_structures/dynamic_array.hpp"

namespace parkway {
namespace data_structures {

/* requires positive int keys */
class map_to_pos_int {
 protected:
  int numEntries;
  int size;
  int useHash;

  DynamicArray<int> entries;
  DynamicArray<int> table;
  DynamicArray<int> keys;

 public:
  map_to_pos_int();
  map_to_pos_int(int _size, int use_hash);
  ~map_to_pos_int() {}

  void createTable(int _size, int use_hash);
  void destroyTable();
  void recoverTable();

  int insertKey(int key, int val);
  void resetSlots();
  int insertIfEmpty(int key, int val);

  int getCareful(int key);
  int getVal(int key);

  inline int usingHashTable() {
    return useHash;
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

#endif  // _DATA_STRUCTURES_MAP_TO_POS_INT_HPP
