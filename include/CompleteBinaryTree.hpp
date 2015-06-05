
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

using namespace std;

template <class T> class CompleteBinaryTree {
protected:
  FastDynaArray<int> tree;
  FastDynaArray<T> roots;

  int numRoots;

public:
  CompleteBinaryTree();
  CompleteBinaryTree(const T *rootVals, const int *cmprs, int nRoots);
  ~CompleteBinaryTree() {}

  inline int getRootVal(register int vID) const {
    register int i = 1;

    while (i < numRoots) {
      if (vID < tree[i])
        i = Shiftl(i, 1);
      else
        i = Or(Shiftl(i, 1), 0x1);
    }

    return roots[i - numRoots];
  }

  void setupTree(const T *rootVals, const int *cmprs, int nRoots);
  void fill(const int *cmprs);
};

#endif
