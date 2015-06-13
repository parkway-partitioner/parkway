
#ifndef _PARA_GREEDYKWAY_REFINER_HPP
#define _PARA_GREEDYKWAY_REFINER_HPP

// ### ParaGreedyKwayRefiner.hpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 30/11/2004: Last Modified
//
// ###

#include <iostream>
#include "data_structures/bit_field.hpp"
#include "data_structures/movement_set_table.hpp"
#include "hypergraph/parallel/hypergraph.hpp"
#include "ParaRefiner.hpp"

using parkway::hypergraph::parallel::hypergraph;
namespace ds = parkway::data_structures;

class ParaGreedyKwayRefiner : public ParaRefiner {
protected:
  int numTotVerticesMoved;
  int earlyExit;

  double limit;

  // data_ structures from point of view of vertices

  ds::dynamic_array<int> numNeighParts;
  ds::dynamic_array<int> neighboursOfV;
  ds::dynamic_array<int> neighboursOfVOffsets;

  // data_ structures from point of view of hyperedges

  ds::dynamic_array<int> hEdgeVinPart;
  ds::dynamic_array<int> hEdgeVinPartOffsets;

  // auxiliary structures

  ds::dynamic_array<int> vertices;
  ds::dynamic_array<int> movedVertices;
  ds::dynamic_array<int> seenVertices;
  ds::dynamic_array<int> numPartsSpanned;
  ds::dynamic_array<int> spannedParts;

  bit_field locked;
  bit_field vertSeen;

  // move set structures

  ds::dynamic_array<ds::dynamic_array<int> *> moveSets;
  ds::dynamic_array<int> moveSetData;
  ds::dynamic_array<int> indexIntoMoveSetData;
  ds::dynamic_array<int> numVerticesMoved;

  ds::movement_set_table *movementSets;

public:
  ParaGreedyKwayRefiner(int rank, int nProcs, int nParts, int numVperP,
                        int eExit, double lim, std::ostream &out);
  ~ParaGreedyKwayRefiner();

  void dispRefinementOptions() const;
  void release_memory();
  void initDataStructs(const hypergraph &h, MPI_Comm comm);
  void resetDataStructs();
  void setPartitioningStructs(int pNumber, MPI_Comm comm);
  // void initVertexPartTable(MPI_Comm comm);
  void refine(hypergraph &h, MPI_Comm comm);

  int runGreedyKwayRefinement(hypergraph &h, int pNo, MPI_Comm comm);
  int doGreedyPass(int lowToHigh, MPI_Comm comm);
  int computeCutsize(MPI_Comm comm);

  void manageBalanceConstraint(MPI_Comm comm);
  void takeBackPassMoves();

  void updateVertexMoveInfo(MPI_Comm comm);
  void updateAdjVertStatus(int v, int sP, int bestMove);
  void unmakeMoves(int indexIntoMoveSets, int from, int to);

  void sanityHedgeCheck() const;
  void nonLocVertCheck() const;
};

#endif
