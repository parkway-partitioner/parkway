#ifndef _DATA_STRUCTURES_MAP_FROM_POS_INT_HPP
#define _DATA_STRUCTURES_MAP_FROM_POS_INT_HPP

#include "data_structures/dynamic_array.hpp"

namespace parkway {
namespace data_structures {

/* New MapFromPosInt Class */
/* requires positive int keys */
template<typename T> class map_from_pos_int {
 protected:
  int numEntries;
  int size;

  DynamicArray<T> table;
  DynamicArray<int> keys;

 public:
  map_from_pos_int();
  map_from_pos_int(int _size);
  ~map_from_pos_int() {}

  void createTable(int _size);
  void destroyTable();
  void recoverTable();

  int insertKey(int key, T val);
  T &getVal(int key);

  inline int getNumEntries() { return numEntries; }
  inline int getNumSlots() { return size; }
};

}  // namespace data_structures
}  // namespace parkway

#endif  // _DATA_STRUCTURES_MAP_FROM_POS_INT_HPP
