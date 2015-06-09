#ifndef _DYNA_HPP
#define _DYNA_HPP

// ### Dyna.hpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
// Modified from code by William Knottenbelt
//
// HISTORY:
//
// 20/12/2004: Last Modified
//
// ###

#include <iostream>
#include <cstdlib>
#include <cassert>
#include <cstring>
#include <cstdio>
#include "Macros.h"
#include "DynaMem.hpp"

namespace parkway {
namespace data_structures {

template <class T> class DynamicArray {
 public:
  inline DynamicArray() {
    length = 0;
    array = nullptr;
  }

  inline DynamicArray(int size) {
    if (size <= 0) {
      length = 0;
      array = nullptr;
      return;
    }
    array = new T[size];

    if (!array) {
#ifdef DEBUG_BASICS
      char message[512];
      sprintf(message,
              "DynamicArray: (%d): memory allocation (size %d) failed!\n",
              (int)this, size);
      std::cout << message;
#endif
      abort();
    }

    length = size;
    return;
  }

  inline ~DynamicArray() {
    DynaMem<T>::deleteArr(array);
    length = 0;
  }

  inline int getLength() const {
    return length;
  }

  inline int setLength(int size) {
    if (size == length)
      return 1;
    if (!adjust(size - length))
      return 0;
    length = size;
    return 1;
  }

  inline int search(T target) const {
    for (int n = 0; n < length; ++n) {
      if (array[n] == target)
        return n;
    }
    return -1;
  }

  inline T *getArray() const { return array; }

  inline void setArray(T *newArray, int size) {
    DynaMem<T>::deleteArr(array);
    array = newArray;
    length = size;
  }

  inline T &operator[](int index) const {
#ifdef DEBUG_BASICS
    if ((index >= length) || (index < 0)) {
#ifdef DEBUG_BASICS
      char message[512];
      sprintf(message,
              "DynamicArray (%d): Array subscript [%d] out of range (Valid "
              "range is 0 to %d)",
              (int)this, index, length - 1);
      std::cout << message;
      delayT();
#endif
      abort();
    }
#endif

    return array[index];
  }

  inline T &operator()(int index) const {
#ifdef DEBUG_BASICS
    if ((index >= length) || (index < 0)) {
      char message[512];
      sprintf(message,
              "DynamicArray: (%d) Array subscript [%d] out of range (Valid "
              "range is 0 to %d)",
              (int)this, index, length - 1);
      std::cout << message;
      abort();
    }
#endif

    return array[index];
  }

  inline void check(int index) {
    if (index >= length)
      expand(index);
  }

  inline void assign(int index, const T value) {
    if (index >= length)
      expand(index);
    array[index] = value;
  }

  inline int adjust(int size) {
    int newLength = length + size;

    if (newLength <= 0) {
      if (length)
        DynaMem<T>::deleteArr(array);

      length = 0;
      return 1;
    }

    T *extendedArray = new T[size + length];

    if (!extendedArray) {
#ifdef DEBUG_BASICS
      char message[512];
      sprintf(message,
              "DynamicArray: (%d): memory allocation (size %d) failed!\n",
              (int)this, size);
      std::cout << message;
#endif
      abort();
    }

    int i;

    if (size >= 0)
      for (i = 0; i < length; ++i)
        extendedArray[i] = array[i];
    else
      for (i = 0; i < newLength; ++i)
        extendedArray[i] = array[i];

    if (length)
      DynaMem<T>::deleteArr(array);

    length = newLength;
    array = extendedArray;

    return 1;
  }

 protected:
  int length;
  T *array;

  inline int expand(int index) {
    int targetLength = 1;
    while (index >= targetLength)
      targetLength <<= 1;
    if (!adjust(targetLength - length))
      return 0;
    return 1;
  }

  inline void delayT() const {
    for (int i = 0; i < 0xFFFFFFA; i++) {
      int *a = new int;
      delete a;
    }
  }



};

}  // namespace data_structures
}  // namespace parkway

#endif
