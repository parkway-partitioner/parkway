
#ifndef _COMPLETE_BINARY_TREE_HPP
#define _COMPLETE_BINARY_TREE_HPP

// ### CompleteBinaryTree.hpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 1/2/2005: Last Modified
//
// ###

#include "Funct.hpp"
#include "data_structures/DynamicArray.h"

namespace parkway {
namespace data_structures {

template <typename T> class CompleteBinaryTree {
 protected:
  DynamicArray<int> tree;
  DynamicArray<T> roots;

  int numRoots;

 public:
  CompleteBinaryTree() : numRoots(0) {
    tree.setLength(0);
    roots.setLength(0);
  }

  CompleteBinaryTree(const T *rootVals, const int *cmprs, int nRoots) {
    setupTree(rootVals, cmprs, nRoots);
  }

  ~CompleteBinaryTree() {
  }

  inline int getRootVal(int vID) const {
    int i = 1;
    while (i < numRoots) {
      if (vID < tree[i]) {
        i = i << 1;
      } else {
        i = (i << 1) | 0x1;
      }
    }
    return roots[i - numRoots];
  }

  void setupTree(const T *rootVals, const int *cmprs, int num_roots) {
    numRoots = num_roots;
    tree.setLength(numRoots);
    roots.setLength(numRoots);

    for (int i = 0; i < numRoots; ++i) {
      roots[i] = rootVals[i];
    }

    tree[0] = 0;
    tree[1] = cmprs[numRoots >> 1];

    if (Funct::log2(numRoots) > 1) {
      fill(cmprs);
    }
  }

  void fill(const int *cmprs) {
    int search_len = numRoots - 1;
    int current_num;
    int j = 0;
    int ij;
    int jk;
    int half = 0;

    DynamicArray<int> search_list(search_len);

    search_list[j++] = numRoots >> 1;

    for (int i = 0; i < (numRoots >> 1) - 1;) {
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
      tree[j++] = cmprs[search_list[i]];
    }
  }


 };

}  // namespace data_structures
}  // namespace parkway

#endif
