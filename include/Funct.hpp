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
#include "data_structures/dynamic_array.hpp"

#ifdef USE_SPRNG
#define SIMPLE_SPRNG
#define USE_MPI
#include "sprng.h"
#endif

/* useful macros */

#define RotateLeft(a, b)                                                       \
  (b > (sizeof(HashKey)) ? (0) : ( ((a >> (sizeof(HashKey) - b))) | (a << b) ))

template <typename T>
T RANDOM(T a, T b) {
#ifdef USE_SPRNG
  return a == b ? a : (static_cast<int>(sprng() * ((b) - (a))) + (a));
#else
  return a == b ? a : (static_cast<int>(drand48() * ((b) - (a))) + (a));
#endif
}

using parkway::data_structures::dynamic_array;

class Funct {
  static int shift1;
  static int shift2;
  static int maxHedgeLen;

  static int tableSizes[16];

 public:
  Funct();
  ~Funct();

  inline static void setMaxHedgeLen(int m) { maxHedgeLen = m; }

  static void printIntro();
  static void printEnd();

  static int setTableSize(int approxNumElem);
  static int getParameterAsInteger(int argc, char **argv, const char *cmpr,
                                   int def);
  static char *getParameterAsCharPtr(int argc, char **argv, const char *cmpr,
                                     char *def);
  static double getParameterAsDouble(int argc, char **argv, const char *cmpr,
                                     double def);
  static double toRecurBal(double e, int nP);

  static HashKey computeHash(const int *vs, int len);
};

#endif
