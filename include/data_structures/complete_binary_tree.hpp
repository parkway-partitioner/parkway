#ifndef _DATA_STRUCTURES_COMPLETE_BINARY_TREE_HPP
#define _DATA_STRUCTURES_COMPLETE_BINARY_TREE_HPP

// ### complete_binary_tree.hpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 1/2/2005: Last Modified
//
// ###

// TODO(gb610): work out what is going on in here; doesn't seem to work for low
// indices -- is this by design or not?

#include "Funct.hpp"
#include "data_structures/dynamic_array.hpp"
#include "utility/math.hpp"

namespace parkway {
namespace data_structures {

template <typename Type> class complete_binary_tree {
 protected:
  dynamic_array<int> keys_;
  dynamic_array<Type> values_;

  int size_;

 public:
  complete_binary_tree() : keys_(0), values_(0), size_(0) {
  }

  complete_binary_tree(const Type *values, const int *keys, int entries) {
    setup(values, keys, entries);
  }

  ~complete_binary_tree() {
  }

  inline int size() const {
    return size_;
  }

  inline int root_value(int vertex_id) const {
    int i = 1;
    while (i < size_) {
      if (vertex_id < keys_[i]) {
        i = i << 1;
      } else {
        i = (i << 1) | 0x1;
      }
    }
    return values_[i - size_];
  }

  void setup(const Type *values, const int *keys, int entries) {
    size_ = entries;
    keys_.reserve(size_);
    values_.reserve(size_);

    for (int i = 0; i < size_; ++i) {
      values_[i] = values[i];
    }

    keys_[0] = 0;
    keys_[1] = keys[size_ >> 1];

    if (utility::math::log2(size_) > 1) {
      fill(keys);
    }
  }

  void fill(const int *keys) {
    int search_len = size_ - 1;
    int current_num;
    int j = 0;
    int ij;
    int jk;
    int half = 0;

    dynamic_array<int> search_list(search_len);

    search_list[j++] = size_ >> 1;

    for (int i = 0; i < (size_ >> 1) - 1;) {
      current_num = search_list[i++];

      if (utility::math::is_power_of_2(current_num)) {
        half = current_num >> 1;
        ij = half;
        jk = current_num + half;
      } else {
        ij = current_num - half;
        jk = current_num + half;
      }

      search_list[j++] = ij;
      search_list[j++] = jk;
    }

    j = 1;
    for (int i = 0; i < search_len; ++i) {
      keys_[j++] = keys[search_list[i]];
    }
  }
 };

}  // namespace data_structures
}  // namespace parkway

#endif
