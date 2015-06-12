
#ifndef _KWAY_GREEDY_REFINER_HPP
#define _KWAY_GREEDY_REFINER_HPP

// ### GreedyKwayRefiner.hpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 25/4/2004: Last Modified
//
// ###

#include "Refiner.hpp"
#include "hypergraph/VertexNode.hpp"

using namespace std;

class GreedyKwayRefiner : public Refiner {

protected:
  int numNonPosMoves;
  double limit;

  // ###
  // data_ structures from point of view of vertices
  // ###

  dynamic_array<int> numNeighParts;
  dynamic_array<int> neighboursOfV;
  dynamic_array<int> neighboursOfVOffsets;

  // ###
  // data_ structures from point of view of hyperedges
  // ###

  dynamic_array<int> hEdgeVinPart;
  dynamic_array<int> hEdgeVinPartOffsets;

  // ###
  // auxiliary structures
  // ###

  dynamic_array<int> vertices;
  dynamic_array<int> vertSeen;
  dynamic_array<int> seenVertices;
  dynamic_array<int> partsSpanned;

public:
  GreedyKwayRefiner(int max, int nparts, double ave, double limit, int dL);
  ~GreedyKwayRefiner();

  void dispRefinerOptions(ostream &out) const;
  void buildDataStructs();
  void destroyDataStructs();
  int initDataStructs();
  void updateAdjVertStat(int v, int sP, int bestDP);

  void refine(serial_hypergraph &h);
  void rebalance(serial_hypergraph &h);

  int runGreedyPass();
  int runRebalancingPass();

  inline void setLimit() {
    double a = static_cast<double>(numVertices) * limit;
    numNonPosMoves = static_cast<int>(floor(a));
  }

  inline int findHeaviestOverweight() const {
    int i = 0;
    int p = -1;

    for (; i < numParts; ++i) {
      if (partWeights[i] > maxPartWt) {
        if (p == -1)
          p = i;
        else if (partWeights[i] > partWeights[p])
          p = i;
      }
    }

    return p;
  }
};

#endif
