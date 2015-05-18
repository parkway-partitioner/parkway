
#  ifndef _DYNA_MEM_HPP
#  define _DYNA_MEM_HPP


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


using namespace std;


template <class T>
class DynaMem
{

public:

  static inline void deletePtr(T* &ptr)
  {
    if(ptr)
      delete ptr;
    ptr = NULL;
  }
  
  static inline void deleteArr(T* &ptr)
  {
    if(ptr)
      delete [] ptr;
    ptr = NULL;
  }

};



#  endif
