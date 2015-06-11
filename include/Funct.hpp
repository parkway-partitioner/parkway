#ifndef _FUNCT_HPP
#define _FUNCT_HPP

// ### Funct.hpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 12/1/2005: Last Modified
//
// NOTES:
//
// - shift1 and shift2 need to
//   be relatively prime to sizeof(HashKey)
//
// ###

#include <iostream>
#include <cmath>
#include <fstream>
#include <unistd.h>
#include <cstdio>
#include <cstdlib>

#include "Macros.h"
#include "Log.h"
#include "data_structures/dynamic_array.hpp"

#ifdef USE_SPRNG
#define SIMPLE_SPRNG
#define USE_MPI
#include "sprng.h"
#endif

/* useful macros */

#define RotateLeft(a, b)                                                       \
  ((b) > (sizeof(HashKey)) ? (0) : (Or(Shiftr((a), sizeof(HashKey) - (b)),     \
                                       Shiftl((a), (b)))))
#ifdef USE_SPRNG
#define RANDOM(a, b)                                                           \
  ((a) == (b) ? (a) : (static_cast<int>((sprng()) * ((b) - (a))) + (a)))
#else
#define RANDOM(a, b)                                                           \
  ((a) == (b) ? (a) : (static_cast<int>(drand48() * ((b) - (a))) + (a)))
#endif

using parkway::data_structures::dynamic_array;

class Funct {
  static char *startHeap;
  static int shift1;
  static int shift2;
  static int maxHedgeLen;

  static int tableSizes[16];

 public:
  Funct();
  ~Funct();

  inline static int getShift1() { return shift1; }
  inline static int getShift2() { return shift2; }
  inline static int getMaxHedgeLen() { return maxHedgeLen; }

  inline static void setShift1(int s) { shift1 = s; }
  inline static void setShift2(int s) { shift2 = s; }
  inline static void setMaxHedgeLen(int m) { maxHedgeLen = m; }

  inline static int search(const int *array, const int length,
                           const int target) {
    int i = 0;

    for (; i < length; ++i)
      if (array[i] == target)
        return i;

    return -1;
  }

  inline static void qsort(const int left, const int right, int *array) {
    int left_arrow = left;
    int right_arrow = right;

    int pivot = array[(left + right) / 2];

    do {
      while (array[right_arrow] > pivot) {
        --right_arrow;
      }
      while (array[left_arrow] < pivot) {
        ++left_arrow;
      }

      if (left_arrow <= right_arrow) {
        std::swap(array[left_arrow], array[right_arrow]);
        ++left_arrow;
        --right_arrow;
      }
    } while (right_arrow >= left_arrow);

    if (left < right_arrow) {
      qsort(left, right_arrow, array);
    }

    if (left_arrow < right) {
      qsort(left_arrow, right, array);
    }
  }

  inline static int log2(int a) {
    int result = 0;
    while (a >>= 1) {
      ++result;
    }
    return result;
  }

  inline static int isPowerOf2(int num) {
    if (num == 1)
      return 1;

    while (num > 1) {
      if (And(num, 0x1))
        return 0;
      num = Shiftr(num, 1);
    }

    return 1;
  }

  static void printIntro(std::ostream &out);
  static void printEnd(std::ostream &out);
  static void printMemUse(int myRank, const char *info);

  static void initStartMem();
  static void write_to_memory();
  static void randomPermutation(int *array, int size);
  static void qsortByAnotherArray(const int left, const int right, int *array,
                                  const int *valArray, const int order);

  static int setTableSize(int approxNumElem);
  static int getParameterAsInteger(int argc, char **argv, const char *cmpr,
                                   int def);
  static char *getParameterAsCharPtr(int argc, char **argv, const char *cmpr,
                                     char *def);
  static double getParameterAsDouble(int argc, char **argv, const char *cmpr,
                                     double def);
  static double toRecurBal(double e, int nP);

  static HashKey computeHash(const int *vs, int len);
  static int getTableSize(int n);
};

#endif
