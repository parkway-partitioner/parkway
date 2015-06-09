
#ifndef _COMPLETE_BINARY_TREE_HPP
#define _COMPLETE_BINARY_TREE_HPP

// ### complete_binary_tree.hpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 1/2/2005: Last Modified
//
// ###

#include "Funct.hpp"
#include "data_structures/dynamic_array.hpp"

namespace parkway {
namespace data_structures {

template <typename T> class complete_binary_tree {
 protected:
  dynamic_array<int> tree_;
  dynamic_array<T> roots_;

  int number_of_roots_;

 public:
  complete_binary_tree() : number_of_roots_(0) {
    tree_.reserve(0);
    roots_.reserve(0);
  }

  complete_binary_tree(const T *rootVals, const int *cmprs, int nRoots) {
    setup(rootVals, cmprs, nRoots);
  }

  ~complete_binary_tree() {
  }

  inline int root_value(int vID) const {
    int i = 1;
    while (i < number_of_roots_) {
      if (vID < tree_[i]) {
        i = i << 1;
      } else {
        i = (i << 1) | 0x1;
      }
    }
    return roots_[i - number_of_roots_];
  }

  void setup(const T *rootVals, const int *cmprs, int num_roots) {
    number_of_roots_ = num_roots;
    tree_.reserve(number_of_roots_);
    roots_.reserve(number_of_roots_);

    for (int i = 0; i < number_of_roots_; ++i) {
      roots_[i] = rootVals[i];
    }

    tree_[0] = 0;
    tree_[1] = cmprs[number_of_roots_ >> 1];

    if (Funct::log2(number_of_roots_) > 1) {
      fill(cmprs);
    }
  }

  void fill(const int *cmprs) {
    int search_len = number_of_roots_ - 1;
    int current_num;
    int j = 0;
    int ij;
    int jk;
    int half = 0;

    dynamic_array<int> search_list(search_len);

    search_list[j++] = number_of_roots_ >> 1;

    for (int i = 0; i < (number_of_roots_ >> 1) - 1;) {
      current_num = search_list[i++];

      if (Funct::isPowerOf2(current_num)) {
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
      tree_[j++] = cmprs[search_list[i]];
    }
  }


 };

}  // namespace data_structures
}  // namespace parkway

#endif
