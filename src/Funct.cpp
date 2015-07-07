//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 12/1/2005: Last Modified
//
// ###
#include "Funct.hpp"
#include "utility/logging.hpp"
#include "utility/math.hpp"
#include "configuration.hpp"
#include <cstring>
#include <cassert>

int Funct::shift1 = 7;
int Funct::shift2 = 13;
int Funct::maxHedgeLen = 0;

Funct::Funct() {}

Funct::~Funct() {}

int Funct::setTableSize(int approxNumElem) {
  if (approxNumElem < VSMALL)
    return VSMALL_TABLE;

  if (approxNumElem < SMALL)
    return SMALL_TABLE;

  if (approxNumElem < MEDIUM)
    return MEDIUM_TABLE;

  if (approxNumElem < LARGE)
    return LARGE_TABLE;

  if (approxNumElem < VLARGE)
    return VLARGE_TABLE;

  else
    return HUGE_TABLE;
}

int Funct::getParameterAsInteger(int argc, char **argv, const char *cmpr,
                                 int def) {
  int i;
  int j = argc - 1;

  for (i = 0; i < j; ++i) {
    if (strcmp(cmpr, argv[i]) == 0)
      return atoi(argv[i + 1]);
  }

  return def;
}

double Funct::getParameterAsDouble(int argc, char **argv, const char *cmpr,
                                   double def) {
  int i;
  int j = argc - 1;

  for (i = 0; i < j; ++i) {
    if (strcmp(cmpr, argv[i]) == 0)
      return atof(argv[i + 1]);
  }

  return def;
}

char *Funct::getParameterAsCharPtr(int argc, char **argv, const char *cmpr,
                                   char *def) {
  int i;
  int j = argc - 1;

  for (i = 0; i < j; ++i) {
    if (strcmp(cmpr, argv[i]) == 0)
      return argv[i + 1];
  }

  return def;
}

double Funct::toRecurBal(double e, int nP) {
  double a;
  double b;

  if (parkway::utility::math::is_power_of_2(nP)) {
    a = (1.0 + e) / nP;
    b = 1.0 / parkway::utility::math::log2(nP);
  } else {
    /* give the balance for nearest
       power of 2 above nP          */

    int pow2 = 1;
    int index = 1;

    while (index < nP) {
      index = index << 1;
      ++pow2;
    }

    a = (1.0 + e) / index;
    b = 1.0 / pow2;
  }

  return (pow(a, b) - 0.5);
}

HashKey Funct::computeHash(const int *vs, int len) {
  HashKey key = 0;

  unsigned int slide1 = SLIDE1;
  unsigned int slide2 = SLIDE2;

  int sum = 0;
  int i;

  for (i = 0; i < maxHedgeLen; ++i) {
    if (i < len) {
      sum += vs[i];
      key ^= RotateLeft(vs[i], slide1);
    } else {
      sum += 1;
      key ^= RotateLeft(1, slide1);
    }

    key ^= RotateLeft(sum, slide2);

    slide1 = (slide1 + shift1) % sizeof(HashKey);
    slide2 = (slide2 + shift2) % sizeof(HashKey);
  }

  return key;
}


void Funct::printIntro() {
  info("\n ------- Parkway v%s -------\n|\n", VERSION_STRING);
}

void Funct::printEnd() {
  info(" --------------------------\n\n");
}
