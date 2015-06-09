
#ifndef _LOADER_HPP
#define _LOADER_HPP

// ### HypergraphLoader.hpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 4/1/2005: Last Modified
//
// ###

#include "Hypergraph.hpp"
#include "data_structures/bit_field.hpp"

using namespace parkway::data_structures;

class HypergraphLoader {
protected:
  int dispOption;
  int currPercentile;

  int numHedges;
  int numVertices;
  int numPins;
  int numPartitions;

  int *vWeight;
  int *hEdgeWeight;
  int *matchVector;
  int *pinList;
  int *hEdgeOffsets;
  int *vToHedges;
  int *vOffsets;

  int *partitionVectors;
  int *partitionOffsets;
  int *partitionCutsizes;

  inline void loadGraph(const Hypergraph &h) {
    numVertices = h.getNumVertices();
    numHedges = h.getNumHedges();
    numPins = h.getNumPins();

    vWeight = h.getVerWeightsArray();
    hEdgeWeight = h.getHedgeWeightsArray();
    pinList = h.getPinListArray();
    hEdgeOffsets = h.getHedgeOffsetArray();
    vToHedges = h.getVtoHedgesArray();
    vOffsets = h.getVerOffsetsArray();
  }

public:
  HypergraphLoader(int disp);
  ~HypergraphLoader();

  void computeHedgesToLoad(bit_field &toLoad);

  inline void loadHypergraphForCoarsening(const Hypergraph &h) {
    loadGraph(h);
    matchVector = h.getMatchVectorArray();
  }

  inline void loadHypergraphForRestrCoarsening(const Hypergraph &h) {
    loadGraph(h);

    matchVector = h.getMatchVectorArray();
    numPartitions = h.getNumPartitions();
    partitionVectors = h.getPartVectorArray();
    partitionOffsets = h.getPartOffsetArray();
    partitionCutsizes = h.getPartCutArray();
  }

  inline void loadHypergraphForRefinement(const Hypergraph &h) {
    loadGraph(h);

    numPartitions = h.getNumPartitions();
    partitionVectors = h.getPartVectorArray();
    partitionOffsets = h.getPartOffsetArray();
    partitionCutsizes = h.getPartCutArray();
  }

  inline void loadHypergraphForSplitting(const Hypergraph &h) {
    loadGraph(h);

    numPartitions = 1;
    partitionVectors = h.getPartVectorArray();
    partitionOffsets = h.getPartOffsetArray();
    partitionCutsizes = h.getPartCutArray();
  }

  inline int getPercentile() const { return currPercentile; }
  inline void setPercentile(int p) { currPercentile = p; }
};

#endif
