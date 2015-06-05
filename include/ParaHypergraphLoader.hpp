
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
#include "Bit.hpp"

using namespace std;

class ParaHypergraphLoader : public GlobalCommunicator {

protected:
  ostream &out_stream;

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

  FastDynaArray<int> hEdgeWeight;
  FastDynaArray<int> hEdgeOffset;
  FastDynaArray<int> locPinList;

  FastDynaArray<int> vToHedgesOffset;
  FastDynaArray<int> vToHedgesList;

  FastDynaArray<int> allocHedges;

public:
  ParaHypergraphLoader(int rank, int nProcs, int nParts, ostream &o);

  virtual ~ParaHypergraphLoader();
  virtual void releaseMemory() = 0;
  virtual void loadHyperGraph(const ParaHypergraph &h, MPI_Comm comm) = 0;

  void computeHedgesToLoad(BitField &toLoad, int numH, int numLocalPins,
                           int *hEdgeWts, int *hEdgeOffsets, MPI_Comm comm);

  inline int getNumParts() const { return numParts; }
  inline int getPercentile() const { return currPercentile; }

  inline void setDispOption(register int d) { dispOption = d; }
  inline void setPercentile(register int p) { currPercentile = p; }
};

#endif
