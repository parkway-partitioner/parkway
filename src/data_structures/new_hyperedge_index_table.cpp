#include "data_structures/new_hyperedge_index_table.hpp"
#include "data_structures/internal/table_utils.hpp"

namespace parkway {
namespace data_structures {

new_hyperedge_index_table::new_hyperedge_index_table(unsigned int _size) {
  createTable(_size);
}

void new_hyperedge_index_table::createTable(unsigned int _size) {
  numEntries = 0;
  size = internal::table_utils::table_size(_size);

#ifdef DEBUG_TABLES
  assert(size >= _size);
#endif

  keys.reserve(size);
  table.reserve(size);
  nextSameKey.reserve(size);

  int i;

  for (i = 0; i < size; ++i) {
    table[i] = -1;
    nextSameKey[i] = -1;
  }
}

void new_hyperedge_index_table::destroyTable() {
  numEntries = 0;
  table.reserve(0);
  keys.reserve(0);
  nextSameKey.reserve(0);
}

void new_hyperedge_index_table::recoverTable() {
  table.reserve(size);
  keys.reserve(size);
  nextSameKey.reserve(size);
  numEntries = 0;

  int i;

  for (i = 0; i < size; ++i) {
    table[i] = -1;
    nextSameKey[i] = -1;
  }
}

void new_hyperedge_index_table::insertKey(HashKey key, int index) {
  int slot = internal::hashes::primary(key, size);
  int lastSeen = -1;

  while (table[slot] != -1) {
    if (keys[slot] == key)
      lastSeen = slot;

    slot = internal::hashes::chained(slot, key, size);
  }

  table[slot] = index;
  keys[slot] = key;

  if (lastSeen >= 0)
    nextSameKey[lastSeen] = slot;

  ++numEntries;
#ifdef DEBUG_TABLES
  assert(numEntries < size);
#endif
}

int new_hyperedge_index_table::getHedgeIndex(HashKey key, int &numSeen) {
#ifdef DEBUG_TABLES
  assert(numSeen < size);
#endif

  int slot;

  if (numSeen == -1) {
    slot = internal::hashes::primary(key, size);
    while (table[slot] != -1) {
      if (keys[slot] == key) {
        numSeen = nextSameKey[slot];
        #ifdef DEBUG_TABLES
        assert(numSeen < size);
        #endif
        return table[slot];
      }
      slot = internal::hashes::chained(slot, key, size);
    }
    return -1;
  }

  slot = numSeen;
  numSeen = nextSameKey[slot];
  #ifdef DEBUG_TABLES
  assert(numSeen < size);
  #endif
  return table[slot];
}

void new_hyperedge_index_table::print() {
  for (int i = 0; i < size; ++i) {
    if (table[i] != -1)
      std::cout << "[" << i << "]: "
           << "key = " << keys[i] << " index = " << table[i]
           << " next = " << nextSameKey[i] << std::endl;
    else
      std::cout << "[empty]" << std::endl;
  }
}

}  // namespace data_structures
}  // namespace parkway
