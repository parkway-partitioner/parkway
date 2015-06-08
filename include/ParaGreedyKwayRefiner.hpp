
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

#include "Bit.hpp"
#include "MovementSets.hpp"
#include "ParaRefiner.hpp"

using namespace std;

class ParaGreedyKwayRefiner : public ParaRefiner {
protected:
  int numTotVerticesMoved;
  int earlyExit;

  double limit;

  // data structures from point of view of vertices

  DynamicArray<int> numNeighParts;
  DynamicArray<int> neighboursOfV;
  DynamicArray<int> neighboursOfVOffsets;

  // data structures from point of view of hyperedges

  DynamicArray<int> hEdgeVinPart;
  DynamicArray<int> hEdgeVinPartOffsets;

  // auxiliary structures

  DynamicArray<int> vertices;
  DynamicArray<int> movedVertices;
  DynamicArray<int> seenVertices;
  DynamicArray<int> numPartsSpanned;
  DynamicArray<int> spannedParts;

  BitField locked;
  BitField vertSeen;

  // move set structures

  DynamicArray<DynamicArray<int> *> moveSets;
  DynamicArray<int> moveSetData;
  DynamicArray<int> indexIntoMoveSetData;
  DynamicArray<int> numVerticesMoved;

  MovementSetTable *movementSets;

public:
  ParaGreedyKwayRefiner(int rank, int nProcs, int nParts, int numVperP,
                        int eExit, double lim, ostream &out);
  ~ParaGreedyKwayRefiner();

  void dispRefinementOptions() const;
  void releaseMemory();
  void initDataStructs(const ParaHypergraph &h, MPI_Comm comm);
  void resetDataStructs();
  void setPartitioningStructs(int pNumber, MPI_Comm comm);
  // void initVertexPartTable(MPI_Comm comm);
  void refine(ParaHypergraph &h, MPI_Comm comm);

  int runGreedyKwayRefinement(ParaHypergraph &h, int pNo, MPI_Comm comm);
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
