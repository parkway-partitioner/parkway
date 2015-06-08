
#ifndef _COMPLETE_BINARY_TREE_CPP
#define _COMPLETE_BINARY_TREE_CPP

// ### CompleteBinaryTree.cpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 1/2/2005: Last Modified
//
// ###

#include "CompleteBinaryTree.hpp"

template <class T> CompleteBinaryTree<T>::CompleteBinaryTree() {
  numRoots = 0;

  tree.setLength(0);
  roots.setLength(0);
}

template <class T>
CompleteBinaryTree<T>::CompleteBinaryTree(const T *rootVals, const int *cmprs,
                                          int nRoots) {
  setupTree(rootVals, cmprs, nRoots);
}

template <class T>
void CompleteBinaryTree<T>::setupTree(const T *rootVals, const int *cmprs,
                                      int nRoots) {
  numRoots = nRoots;

  tree.setLength(nRoots);
  roots.setLength(nRoots);

  for (int i = 0; i < numRoots; ++i)
    roots[i] = rootVals[i];

  tree[0] = 0;
  tree[1] = cmprs[Shiftr(numRoots, 1)];

  if (Funct::log2(numRoots) > 1)
    fill(cmprs);
}

template <class T> void CompleteBinaryTree<T>::fill(const int *cmprs) {
  int i = 0;
  int j = 0;
  int searchLen = numRoots - 1;
  int current_num;
  int ij;
  int jk;
  int half = 0;

  DynamicArray<int> searchList(searchLen);

  searchList[j++] = Shiftr(numRoots, 1);

  for (; i < Shiftr(numRoots, 1) - 1;) {
    current_num = searchList[i++];

    if (Funct::isPowerOf2(current_num)) {
      half = Shiftr(current_num, 1);
      ij = half;
      jk = current_num + half;
    } else {
      ij = current_num - half;
      jk = current_num + half;
    }

    searchList[j++] = ij;
    searchList[j++] = jk;
  }

  j = 1;
  for (i = 0; i < searchLen; ++i)
    tree[j++] = cmprs[searchList[i]];
}

template class CompleteBinaryTree<int>;

#endif
