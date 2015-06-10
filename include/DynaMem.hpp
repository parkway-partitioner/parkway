#ifndef _DYNA_MEM_HPP
#define _DYNA_MEM_HPP

// ### DynaMem.hpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
// Modified from code by William Knottenbelt
//
// HISTORY:
//
// 20/12/2004: Last Modified
//
// ###

// for debug purposes, remove dependence on MemoryTracker
//#  include "MemoryTracker.hpp"
#include <memory>

namespace DynaMem {

template <typename T> inline void deletePtr(T *&ptr) {
  if (ptr != nullptr) {
    delete ptr;
  }
  ptr = nullptr;
}

template <typename T> inline void deleteArr(T *&ptr) {
  if (ptr != nullptr) {
    delete[] ptr;
  }
  ptr = nullptr;
}

template <typename T> inline std::shared_ptr<T> shared_array(std::size_t size) {
  return std::shared_ptr<T>(new T[size], std::default_delete<T[]>());
}

}
#endif
