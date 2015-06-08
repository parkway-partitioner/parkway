
#ifndef COUNT_HPP
#define COUNT_HPP

// ### Count.hpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 31/1/2005: Last Modified
//
// ###

class Count {
private:
  static int numNew;
  static int numDel;

public:
  Count() {}
  ~Count() {}

  static inline int getNumDel() { return numDel; }
  static inline int getNumNew() { return numNew; }

  static inline void incDel() { ++numDel; }
  static inline void incNew() { ++numNew; }
};

#endif
