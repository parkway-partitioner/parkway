#ifndef _STACK_HPP
#define _STACK_HPP

// ### stack.hpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 30/11/2004: Last Modified
//
// ###

#include "data_structures/dynamic_array.hpp"

namespace parkway {
namespace data_structures {

template <class T> class stack {
 public:
  inline stack() : data_(0), size_(0) {
  }

  inline ~stack() {
  }

  inline void push(const T elem) {
    data_.assign(size_++, elem);
  }

  inline int size() const {
    return size_;
  }

  inline T pop() {
#ifdef DEBUG_BASICS
    assert(numElem > 0);
#endif
    return data_[--size_];
  }

  inline T top() const {
#ifdef DEBUG_BASICS
    assert(numElem > 0);
#endif
    return data_[size_ - 1];
  }

 protected:
  dynamic_array<T> data_;
  int size_;
};

}  // namespace data_structures
}  // namespace parkway

#endif
