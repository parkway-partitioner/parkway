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
#include <memory>

namespace dynamic_memory {

template <typename T> inline void delete_pointer(T *&ptr) {
  if (ptr != nullptr) {
    delete ptr;
  }
  ptr = nullptr;
}

template <typename T> inline void delete_array(T *&ptr) {
  if (ptr != nullptr) {
    delete[] ptr;
  }
  ptr = nullptr;
}

}
#endif
