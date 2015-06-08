
#ifndef _PRINT_CHARS_HPP
#define _PRINT_CHARS_HPP

// ### PrintHypergraphChars.hpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 5/1/2005: Last Modified
//
// ###

#include <iostream>
#include <fstream>
#include "DynamicArray.h"
#include "StringUtils.hpp"

using namespace std;

class HypergraphChars {
protected:
  int zerosOnDiag;
  int removeSingletons;
  int numVertices;
  int numHedges;
  int numPins;

  char *graphName;

  FastDynaArray<int> vWeights;
  FastDynaArray<int> hEdgeWeights;
  FastDynaArray<int> hEdgeLens;

public:
  HypergraphChars(int remove, int nZ);
  ~HypergraphChars();

  void initStructs(const char *filename);
  void outputChars(ostream &out);
  void qsort(const int l, const int r, int *array, const int *cmp);

  inline void swap(int &a, int &b) {
    int temp = a;
    a = b;
    b = temp;
  }
};

#endif
