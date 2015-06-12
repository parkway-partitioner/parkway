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

#include "hypergraph.hpp"
#include "data_structures/bit_field.hpp"

using namespace parkway::data_structures;

class hypergraph_loader {
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

  inline void load(const hypergraph &h) {
    numVertices = h.number_of_vertices();
    numHedges = h.number_of_hyperedges();
    numPins = h.number_of_pins();

    vWeight = h.vertex_weights();
    hEdgeWeight = h.hyperedge_weights();
    pinList = h.pin_list();
    hEdgeOffsets = h.hyperedge_offsets();
    vToHedges = h.vertex_to_hyperedges();
    vOffsets = h.vertex_offsets();
  }

public:
  hypergraph_loader(int disp);
  ~hypergraph_loader();

  void compute_hyperedges_to_load(bit_field &toLoad);

  inline void load_for_coarsening(const hypergraph &h) {
    load(h);
    matchVector = h.match_vector();
  }

  inline void loadHypergraphForRestrCoarsening(const hypergraph &h) {
    load(h);

    matchVector = h.match_vector();
    numPartitions = h.number_of_partitions();
    partitionVectors = h.partition_vector();
    partitionOffsets = h.partition_offsets();
    partitionCutsizes = h.partition_cuts();
  }

  inline void load_for_refinement(const hypergraph &h) {
    load(h);

    numPartitions = h.number_of_partitions();
    partitionVectors = h.partition_vector();
    partitionOffsets = h.partition_offsets();
    partitionCutsizes = h.partition_cuts();
  }

  inline void load_for_splitting(const hypergraph &h) {
    load(h);

    numPartitions = 1;
    partitionVectors = h.partition_vector();
    partitionOffsets = h.partition_offsets();
    partitionCutsizes = h.partition_cuts();
  }

  inline int get_percentile() const { return currPercentile; }
  inline void set_percentile(int p) { currPercentile = p; }
};

#endif