

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
#include "DynamicArray.h"

template <class T> class Stack {

protected:
  int numElem;
  DynamicArray<T> array;

public:
  inline Stack() {
    array.setLength(0);
    numElem = 0;
  }

  inline ~Stack() {}

  inline void push(const T elem) { array.assign(numElem++, elem); }
  inline int getNumElem() const { return numElem; }

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
};

#endif
