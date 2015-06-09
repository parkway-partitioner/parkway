
#ifndef _BITFIELD_HPP
#define _BITFIELD_HPP

// ### bit_field.hpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
// Modified from code by William Knottenbelt
//
// HISTORY:
//
// 30/11/2004: Last Modified
//
// NOTES:
//
//  - Assume that sizeof(unsigned int) == 4 bytes
//
// ###

#include <iostream>
#include "data_structures/dynamic_array.hpp"

namespace parkway {
namespace data_structures {

class bit_field {
 public:
  bit_field() {
  }

  bit_field(int bits) {
    setLength(bits);
  }

  ~bit_field() {
  }

  inline void check(int bit) {
    int old = data.capacity();
    int bits = bit + 1;
    int target = (bits >> 5) + ((bits & 31) ? 1 : 0);

    data.check(target);
    length = data.capacity();

    if (length > old) {
      for (int n = old; n < length; ++n)
        data[n] = 0;
    }
    bitLength = (length << 5);
  }

  inline int getNumBits() const {
    return bitLength;
  }

  inline int getLength() const {
    return length;
  }

  inline unsigned int *getData() const {
    return data.data();
  }

  inline unsigned int getChunk(int i) const {
    return data[i];
  }

  inline void setLength(int bits) {
    bitLength = bits;
    length = (bits >> 5) + ((bits & 31) ? 1 : 0);
    data.reserve(length);
  }

  inline void unset() {
    for (int l = 0; l < length; ++l) {
      data[l] = 0;
    }
  }

  inline void set() {
    for (int l = 0; l < length; ++l) {
      data[l] = 0xFFFFFFFF;
    }
  }

  inline bool test_all() {
    bool all_set = true;
    for (int i = 0; i < length; ++i)
      if (data[i] != 0xFFFFFFFF)
        return false;
    return all_set;
  }

  inline int operator()(int index) const {
    return (data(index >> 5) & (1 << (index & 31))) ? 1 : 0;
  }

  inline void set(int index) {
    data[index >> 5] |= (1 << (index & 31));
  }

  inline void unset(int index) {
    data[index >> 5] &= (~(1 << (index & 31)));
  }

 private:
  int bitLength;
  int length;
  dynamic_array<unsigned int> data;
};

}  // namespace data_structures
}  // namespace parkway

#endif
