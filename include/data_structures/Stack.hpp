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
  inline stack() : array(0), numElem(0) {
  }

  inline ~stack() {
  }

  inline void push(const T elem) {
    array.assign(numElem++, elem);
  }

  inline int getNumElem() const {
    return numElem;
  }

  inline T pop() {
#ifdef DEBUG_BASICS
    assert(numElem > 0);
#endif
    return array[--numElem];
  }

  inline T getTopElem() const {
#ifdef DEBUG_BASICS
    assert(numElem > 0);
#endif
    return array[numElem - 1];
  }

 protected:
  dynamic_array<T> array;
  int numElem;
};

}  // namespace data_structures
}  // namespace parkway

#endif
