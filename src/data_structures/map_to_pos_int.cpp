#include "data_structures/map_to_pos_int.hpp"
#include "data_structures/internal/table_utils.hpp"

namespace parkway {
namespace data_structures {

map_to_pos_int::map_to_pos_int() {
  createTable(0, 0);
}

map_to_pos_int::map_to_pos_int(int _size, int use_hash) {
  createTable(_size, use_hash);
}

void map_to_pos_int::createTable(int _size, int use_hash) {
  useHash = use_hash;
  numEntries = 0;
  entries.setLength(2048);

  int i;

  if (useHash == 0) {
    size = _size;
    keys.setLength(0);
    table.setLength(size);

    for (i = 0; i < size; ++i)
      table[i] = -1;
  } else {
    size = internal::table_utils::table_size(_size);
#ifdef DEBUG_TABLES
    assert(size >= _size);
#endif

    keys.setLength(size);
    table.setLength(size);

    for (i = 0; i < size; ++i) {
      keys[i] = -1;
      table[i] = -1;
    }
  }
}

void map_to_pos_int::destroyTable() {
  numEntries = 0;
  entries.setLength(0);
  table.setLength(0);
  keys.setLength(0);
}

void map_to_pos_int::recoverTable() {
  int i;

  numEntries = 0;
  entries.setLength(2048);

  if (useHash) {
    table.setLength(size);
    keys.setLength(size);

    for (i = 0; i < size; ++i) {
      keys[i] = -1;
      table[i] = -1;
    }
  } else {
    table.setLength(size);
    keys.setLength(0);

    for (i = 0; i < size; ++i)
      table[i] = -1;
  }
}

int map_to_pos_int::insertKey(int key, int val) {
  #ifdef DEBUG_TABLES
  assert(key >= 0);
  assert(val >= 0);
  #endif

  if (useHash) {
    int indepKey = internal::table_utils::scatter_key(key);
    int slot = internal::hashes::primary(indepKey, size);

    while (keys[slot] != -1 && keys[slot] != indepKey) {
      slot = (slot + internal::hashes::secondary(indepKey, size)) % size;
    }

    if (keys[slot] == -1) {
      #ifdef DEBUG_TABLES
      assert(table[slot] == -1);
      #endif
      table[slot] = val;
      keys[slot] = indepKey;
      entries.assign(numEntries++, slot);
      return 0;
    } else {
      #ifdef DEBUG_TABLES
      assert(keys[slot] == indepKey);
      assert(table[slot] >= 0);
      #endif
      table[slot] = val;
      return 1;
    }
  } else {
    #ifdef DEBUG_TABLES
    assert(key < size);
    #endif
    if (table[key] == -1) {
      table[key] = val;
      entries.assign(numEntries++, key);
      return 0;
    } else {
      table[key] = val;
      return 1;
    }
  }
}

int map_to_pos_int::insertIfEmpty(int key, int val) {
  #ifdef DEBUG_TABLES
  assert(key >= 0);
  #endif

  if (useHash) {
    int indepKey = internal::table_utils::scatter_key(key);
    int slot = internal::hashes::primary(indepKey, size);

    while (keys[slot] != -1 && keys[slot] != indepKey) {
      slot = (slot + internal::hashes::secondary(indepKey, size)) % size;
    }

    if (keys[slot] == -1) {
      table[slot] = val;
      keys[slot] = indepKey;
      entries.assign(numEntries++, slot);
      return -1;
    } else {
      return table[slot];
    }
  } else {
    #ifdef DEBUG_TABLES
    assert(key < size);
    #endif

    if (table[key] == -1) {
      table[key] = val;
      entries.assign(numEntries++, key);
      return -1;
    } else {
      return table[key];
    }
  }
}

void map_to_pos_int::resetSlots() {
  int i;
  int slot;

  if (useHash) {
    for (i = 0; i < numEntries; ++i) {
      slot = entries[i];

      keys[slot] = -1;
      table[slot] = -1;
    }
  } else {
    for (i = 0; i < numEntries; ++i)
      table[entries[i]] = -1;
  }

  numEntries = 0;
}

int map_to_pos_int::getCareful(int key) {
#ifdef DEBUG_TABLES
  assert(key >= 0);
#endif

  if (useHash) {
    int indepKey = internal::table_utils::scatter_key(key);
    int slot = internal::hashes::primary(indepKey, size);

    while (keys[slot] != indepKey && keys[slot] != -1) {
      slot = (slot + internal::hashes::secondary(indepKey, size)) % size;
    }

    if (keys[slot] == -1) {
      #ifdef DEBUG_TABLES
      assert(table[slot] == -1);
      #endif
      return -1;
    } else {
      #ifdef DEBUG_TABLES
      assert(table[slot] != -1);
      #endif
      return table[slot];
    }
  } else {
    #ifdef DEBUG_TABLES
    assert(key < size);
    #endif
    if (table[key] == -1)
      return -1;
    else
      return table[key];
  }
}

int map_to_pos_int::getVal(int key) {
  #ifdef DEBUG_TABLES
  assert(key >= 0);
  #endif

  if (useHash) {
    int indepKey = internal::table_utils::scatter_key(key);
    int slot = internal::hashes::primary(indepKey, size);

    while (keys[slot] != indepKey) {
      slot = (slot + internal::hashes::secondary(indepKey, size)) % size;
      #ifdef DEBUG_TABLES
      assert(keys[slot] != -1);
      #endif
    }

    return table[slot];
  } else {
    #ifdef DEBUG_TABLES
    assert(key < size);
    #endif
    return table[key];
  }
}

}  // namespace data_structures
}  // namespace parkway
