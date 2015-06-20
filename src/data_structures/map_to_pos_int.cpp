#include "data_structures/map_to_pos_int.hpp"
#include "data_structures/internal/table_utils.hpp"
#include "assert.h"

namespace parkway {
namespace data_structures {

map_to_pos_int::map_to_pos_int() : map_to_pos_int(0, false) {
}

map_to_pos_int::map_to_pos_int(int new_capacity, bool use_hash)
  : use_hash_(use_hash) {
  create(new_capacity, use_hash);
}

void map_to_pos_int::create(int new_capacity, bool use_hash) {
  use_hash_ = use_hash;
  size_ = 0;
  entries.resize(2048);

  if (use_hash_) {
    capacity_ = internal::table_utils::table_size(new_capacity);
    assert(capacity_ >= new_capacity);

    keys.resize(capacity_);
    table.resize(capacity_);

    for (int i = 0; i < capacity_; ++i) {
      keys[i] = -1;
      table[i] = -1;
    }
  } else {
    capacity_ = new_capacity;
    keys.resize(0);
    table.resize(capacity_);

    for (int i = 0; i < capacity_; ++i) {
      table[i] = -1;
    }
  }
}

void map_to_pos_int::destroy() {
  size_ = 0;
  entries.resize(0);
  table.resize(0);
  keys.resize(0);
}

void map_to_pos_int::recover() {
  int i;

  size_ = 0;
  entries.resize(2048);

  if (use_hash_) {
    table.resize(capacity_);
    keys.resize(capacity_);

    for (i = 0; i < capacity_; ++i) {
      keys[i] = -1;
      table[i] = -1;
    }
  } else {
    table.resize(capacity_);
    keys.resize(0);

    for (i = 0; i < capacity_; ++i)
      table[i] = -1;
  }
}

int map_to_pos_int::insert(int key, int val) {
  #ifdef DEBUG_TABLES
  assert(key >= 0);
  assert(val >= 0);
  #endif

  if (use_hash_) {
    int independent_key = internal::table_utils::scatter_key(key);
    int slot = internal::hashes::primary(independent_key, capacity_);

    while (keys[slot] != -1 && keys[slot] != independent_key) {
      slot = (slot + internal::hashes::secondary(independent_key, capacity_)) % capacity_;
    }

    if (keys[slot] == -1) {
      #ifdef DEBUG_TABLES
      assert(table[slot] == -1);
      #endif
      table[slot] = val;
      keys[slot] = independent_key;
      entries[size_++] = slot;
      return false;
    } else {
      #ifdef DEBUG_TABLES
      assert(keys[slot] == independent_key);
      assert(table[slot] >= 0);
      #endif
      table[slot] = val;
      return true;
    }
  } else {
    #ifdef DEBUG_TABLES
    assert(key < capacity_);
    #endif
    if (table[key] == -1) {
      table[key] = val;
      entries[size_++] = key;
      return false;
    } else {
      table[key] = val;
      return true;
    }
  }
}

int map_to_pos_int::insert_if_empty(int key, int val) {
  #ifdef DEBUG_TABLES
  assert(key >= 0);
  #endif

  if (use_hash_) {
    int independent_key = internal::table_utils::scatter_key(key);
    int slot = internal::hashes::primary(independent_key, capacity_);

    while (keys[slot] != -1 && keys[slot] != independent_key) {
      slot = (slot + internal::hashes::secondary(independent_key, capacity_)) % capacity_;
    }

    if (keys[slot] == -1) {
      table[slot] = val;
      keys[slot] = independent_key;
      entries[size_++] = slot;
      return -1;
    } else {
      return table[slot];
    }
  } else {
    #ifdef DEBUG_TABLES
    assert(key < capacity_);
    #endif

    if (table[key] == -1) {
      table[key] = val;
      entries[size_++] = key;
      return -1;
    } else {
      return table[key];
    }
  }
}

void map_to_pos_int::clear() {
  int i;
  int slot;

  if (use_hash_) {
    for (i = 0; i < size_; ++i) {
      slot = entries[i];
      keys[slot] = -1;
      table[slot] = -1;
    }
  } else {
    for (i = 0; i < size_; ++i) {
      table[entries[i]] = -1;
    }
  }

  size_ = 0;
}

int map_to_pos_int::get_careful(int key) {
#ifdef DEBUG_TABLES
  assert(key >= 0);
#endif

  if (use_hash_) {
    int independent_key = internal::table_utils::scatter_key(key);
    int slot = internal::hashes::primary(independent_key, capacity_);

    while (keys[slot] != independent_key && keys[slot] != -1) {
      slot = (slot + internal::hashes::secondary(independent_key, capacity_)) % capacity_;
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
    assert(key < capacity_);
    #endif
    if (table[key] == -1)
      return -1;
    else
      return table[key];
  }
}

int map_to_pos_int::get(int key) {
  #ifdef DEBUG_TABLES
  assert(key >= 0);
  #endif

  if (use_hash_) {
    int independent_key = internal::table_utils::scatter_key(key);
    int slot = internal::hashes::primary(independent_key, capacity_);

    while (keys[slot] != independent_key) {
      slot = (slot + internal::hashes::secondary(independent_key, capacity_)) % capacity_;
      #ifdef DEBUG_TABLES
      assert(keys[slot] != -1);
      #endif
    }

    return table[slot];
  } else {
    #ifdef DEBUG_TABLES
    assert(key < capacity_);
    #endif
    return table[key];
  }
}

}  // namespace data_structures
}  // namespace parkway
