#ifndef _DATA_STRUCTURES_NEW_HYPEREDGE_INDEX_TABLE_HPP
#define _DATA_STRUCTURES_NEW_HYPEREDGE_INDEX_TABLE_HPP

#include "data_structures/dynamic_array.hpp"

namespace parkway {
namespace data_structures {

class new_hyperedge_index_table {
 protected:
  unsigned int numEntries;
  unsigned int size;

  DynamicArray<int> table;
  DynamicArray<int> nextSameKey;
  DynamicArray<HashKey> keys;

 public:
  new_hyperedge_index_table(unsigned int _size);
  ~new_hyperedge_index_table() {}

  void createTable(unsigned int _size);
  void destroyTable();
  void recoverTable();

  void print();
  void insertKey(HashKey key, int index);
  int getHedgeIndex(HashKey key, int &numSeen);

  inline unsigned int getNumEntries() {
    return numEntries;
  }

  inline unsigned int getNumSlots() {
    return size;
  }
};

}  // namespace data_structures
}  // namespace parkway

#endif  // _DATA_STRUCTURES_NEW_HYPEREDGE_INDEX_TABLE_HPP
