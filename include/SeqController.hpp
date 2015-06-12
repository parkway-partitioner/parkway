

#ifndef _SEQ_CONTROLLER_HPP
#define _SEQ_CONTROLLER_HPP

// ### SeqController.hpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 30/11/2004: Last Modified
//
// ###

#include "ParaHypergraph.hpp"
#include "hypergraph.hpp"

using namespace std;

class SeqController {
protected:
  ostream &out_stream;

  int dispOption;
  int myRank;
  int numProcs;
  int numParts;
  int numSeqRuns;
  int maxVertexWt;

  double kWayConstraint;
  double acceptProp;

  hypergraph *h;

  dynamic_array<int> partitionVector;
  dynamic_array<int> partitionCuts;
  dynamic_array<int> partitionVectorOffsets;

public:
  SeqController(int rank, int nProcs, int nParts, ostream &out);

  virtual ~SeqController();
  virtual void dispSeqControllerOptions() const = 0;
  virtual void runSeqPartitioner(ParaHypergraph &hgraph, MPI_Comm comm) = 0;
  virtual void initSeqPartitions(ParaHypergraph &h, MPI_Comm comm);
  virtual void initCoarsestHypergraph(ParaHypergraph &hgraph, MPI_Comm comm);

  int chooseBestPartition() const;
  int getAcceptCut() const;

  inline void setNumSeqRuns(int r) { numSeqRuns = r; }
  inline void setDispOption(int d) { dispOption = d; }
  inline void setMaxVertexWt(int max) { maxVertexWt = max; }
  inline void setHypergraph(hypergraph *hGraph) { h = hGraph; }
  inline void setKwayConstraint(double c) { kWayConstraint = c; }
  inline void setAcceptProp(double p) { acceptProp = p; }
};

#endif
