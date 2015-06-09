#include "data_structures/map_from_pos_int.hpp"
#include "data_structures/internal/table_utils.hpp"

namespace parkway {
namespace data_structures {

template <typename T>
map_from_pos_int<T>::map_from_pos_int() : map_from_pos_int(0) {
}

template <typename T>
map_from_pos_int<T>::map_from_pos_int(int _size) {
  createTable(_size);
}

template <typename T>
void map_from_pos_int<T>::createTable(int _size) {
  numEntries = 0;
  size = internal::table_utils::table_size(_size);

  #ifdef DEBUG_TABLES
  assert(size >= _size);
  #endif

  keys.setLength(size);
  table.setLength(size);

  for (std::size_t i = 0; i < size; ++i) {
    keys[i] = -1;
  }
}

template <typename T>
void map_from_pos_int<T>::destroyTable() {
  numEntries = 0;
  table.setLength(0);
  keys.setLength(0);
}

template <typename T>
void map_from_pos_int<T>::recoverTable() {
  table.setLength(size);
  keys.setLength(size);
  numEntries = 0;

  for (std::size_t i = 0; i < size; ++i) {
    keys[i] = -1;
  }
}

template <typename T>
int map_from_pos_int<T>::insertKey(int key, T val) {
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

template <typename T> T &map_from_pos_int<T>::getVal(int key) {
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


}  // namespace data_structures
}  // namespace parkway
