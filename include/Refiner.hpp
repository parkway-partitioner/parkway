
#ifndef _REFINER_HPP
#define _REFINER_HPP

// ### Refiner.hpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 30/11/2004: Last Modified
//
// ###

#include <ostream>
#include "Macros.h"
#include "HypergraphLoader.hpp"
#include "data_structures/dynamic_array.hpp"

using parkway::data_structures::DynamicArray;

class Refiner : public HypergraphLoader {
 protected:
  int maxPartWt;
  int numParts;
  int *partitionVector;

  double acceptProp;
  double avePartWt;

  DynamicArray<int> partWeights;

 public:
  Refiner(int dispL);

  virtual ~Refiner();
  virtual void refine(Hypergraph &h) = 0;
  virtual void dispRefinerOptions(std::ostream &out) const = 0;

  inline int getMaxPartWt() const { return maxPartWt; }

  inline void setMaxPartWt(int max) { maxPartWt = max; }
  inline void setAvePartWt(double ave) { avePartWt = ave; }
  inline void setPartitionVector(int nP) {
    partitionVector = &partitionVectors[partitionOffsets[nP]];
  }

  int calcCutsize() const;
};

#endif
