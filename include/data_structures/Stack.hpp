

#ifndef _STACK_HPP
#define _STACK_HPP

// ### Stack.hpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 30/11/2004: Last Modified
//
// ###

#include "data_structures/DynamicArray.h"

namespace parkway {
namespace data_structures {

template <class T> class Stack {
 public:
  inline Stack() : array(0), numElem(0) {
  }

  inline ~Stack() {
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
  DynamicArray<T> array;
  int numElem;
};

}  // namespace data_structures
}  // namespace parkway

#endif
