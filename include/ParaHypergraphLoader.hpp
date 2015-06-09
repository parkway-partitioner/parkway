
#ifndef _PARA_HYPO_LOADER_HPP
#define _PARA_HYPO_LOADER_HPP

// ### ParaHypergraphLoader.hpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 4/1/2005: Last Modified
//
// ###

#include "ParaHypergraph.hpp"

#include <iostream>
#include "data_structures/bit_field.hpp"
#include "data_structures/dynamic_array.hpp"

using parkway::data_structures::dynamic_array;
using parkway::data_structures::bit_field;

class ParaHypergraphLoader : public GlobalCommunicator {
 protected:
  std::ostream &out_stream;

  /* hypergraph variables */

  int numParts;
  int numHedges;
  int numLocPins;
  int numLocalVertices;
  int totalVertices;
  int minVertexIndex;
  int maxVertexIndex;
  int locVertWt;

  int numAllocHedges;
  int dispOption;

  /* member for approx coarsening and refinement */

  int currPercentile;

  int *vWeight;
  int *matchVector;

  dynamic_array<int> hEdgeWeight;
  dynamic_array<int> hEdgeOffset;
  dynamic_array<int> locPinList;

  dynamic_array<int> vToHedgesOffset;
  dynamic_array<int> vToHedgesList;

  dynamic_array<int> allocHedges;

 public:
  ParaHypergraphLoader(int rank, int nProcs, int nParts, std::ostream &o);

  virtual ~ParaHypergraphLoader();
  virtual void releaseMemory() = 0;
  virtual void loadHyperGraph(const ParaHypergraph &h, MPI_Comm comm) = 0;

  void computeHedgesToLoad(bit_field &toLoad, int numH, int numLocalPins,
                           int *hEdgeWts, int *hEdgeOffsets, MPI_Comm comm);

  inline int getNumParts() const { return numParts; }
  inline int getPercentile() const { return currPercentile; }

  inline void setDispOption(int d) { dispOption = d; }
  inline void setPercentile(int p) { currPercentile = p; }
};

#endif
