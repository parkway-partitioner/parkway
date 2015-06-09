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

template <class T> class dynamic_array {
 public:
  inline dynamic_array() {
    capacity_ = 0;
    data_ = nullptr;
  }

  inline dynamic_array(int capacity) {
    if (capacity <= 0) {
      capacity_ = 0;
      data_ = nullptr;
      return;
    }
    data_ = new T[capacity];

    if (!data_) {
      abort();
    }

    capacity_ = capacity;
    return;
  }

  inline ~dynamic_array() {
    DynaMem<T>::deleteArr(data_);
    capacity_ = 0;
  }

  inline int capacity() const {
    return capacity_;
  }

  inline bool reserve(int size) {
    if (size == capacity_) {
      return true;
    } else  if (!adjust(size - capacity_)) {
      return false;
    }
    capacity_ = size;
    return true;
  }

  inline int find(T target) const {
    for (int n = 0; n < capacity_; ++n) {
      if (data_[n] == target) {
        return n;
      }
    }
    return -1;
  }

  inline T *data() const {
    return data_;
  }

  inline void set_data(T *new_data, int size) {
    DynaMem<T>::deleteArr(data_);
    data_ = new_data;
    capacity_ = size;
  }

  inline T &operator[](int index) const {
#ifdef DEBUG_BASICS
    if ((index >= capacity_) || (index < 0)) {
#ifdef DEBUG_BASICS
      char message[512];
      sprintf(message,
              "dynamic_array (%d): Array subscript [%d] out of range (Valid "
              "range is 0 to %d)",
              (int)this, index, capacity_ - 1);
      std::cout << message;
      delayT();
#endif
      abort();
    }
#endif

    return data_[index];
  }

  inline T &operator()(int index) const {
#ifdef DEBUG_BASICS
    if ((index >= capacity_) || (index < 0)) {
      char message[512];
      sprintf(message,
              "dynamic_array: (%d) Array subscript [%d] out of range (Valid "
              "range is 0 to %d)",
              (int)this, index, capacity_ - 1);
      std::cout << message;
      abort();
    }
#endif

    return data_[index];
  }

  inline void check(int index) {
    if (index >= capacity_)
      expand(index);
  }

  inline void assign(int index, const T value) {
    if (index >= capacity_)
      expand(index);
    data_[index] = value;
  }

  inline int adjust(int size) {
    int newLength = capacity_ + size;

    if (newLength <= 0) {
      if (capacity_)
        DynaMem<T>::deleteArr(data_);

      capacity_ = 0;
      return 1;
    }

    T *extendedArray = new T[size + capacity_];

    if (!extendedArray) {
#ifdef DEBUG_BASICS
      char message[512];
      sprintf(message,
              "dynamic_array: (%d): memory allocation (size %d) failed!\n",
              (int)this, size);
      std::cout << message;
#endif
      abort();
    }

    int i;

    if (size >= 0)
      for (i = 0; i < capacity_; ++i)
        extendedArray[i] = data_[i];
    else
      for (i = 0; i < newLength; ++i)
        extendedArray[i] = data_[i];

    if (capacity_)
      DynaMem<T>::deleteArr(data_);

    capacity_ = newLength;
    data_ = extendedArray;

    return 1;
  }

 protected:
  int capacity_;
  T *data_;

  inline int expand(int index) {
    int targetLength = 1;
    while (index >= targetLength)
      targetLength <<= 1;
    if (!adjust(targetLength - capacity_))
      return 0;
    return 1;
  }

  inline void delay_time() const {
    for (int i = 0; i < 0xFFFFFFA; i++) {
      int *a = new int;
      delete a;
    }
  }
};

}  // namespace data_structures
}  // namespace parkway

#endif
