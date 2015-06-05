

#ifndef _RECUR_BISECT_CONTROLLER_HPP
#define _RECUR_BISECT_CONTROLLER_HPP

// ### RecurBisectController.hpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 30/11/2004: Last Modified
//
// NOTES:
//
// implement the algorithms such that it is not imperative
// that partition balance is maintained during recursive
// bisection. The greedy k-way refiner can be used to rebalance
// an unbalanced partition
//
// ###

#include "SeqController.hpp"
#include "VCycleFinalBisectionController.hpp"
#include "VCycleAllBisectionController.hpp"
#include "GreedyKwayRefiner.hpp"
#include "Bisection.hpp"

using namespace std;

class RecurBisectController : public SeqController {
protected:
  int numBisectRuns;
  int logK;
  int maxPartWt;
  int sumOfCuts;
  int locVertPartInfoLen;

  int numMyPartitions;

  double bisectConstraint;
  double avePartWt;
  double aveInitBisectionWt;

  BisectionController *bisector;
  GreedyKwayRefiner *kWayRefiner;

  FastDynaArray<int> locVertPartitionInfo;
  FastDynaArray<int> allPartitionInfo;

public:
  RecurBisectController(BisectionController *b, GreedyKwayRefiner *k, int rank,
                        int nProcs, int nParts, int nBisectRuns, ostream &out);
  ~RecurBisectController();

  void dispSeqControllerOptions() const;
  void convToBisectionConstraints();

  void runSeqPartitioner(ParaHypergraph &hgraph, MPI_Comm comm);
  void initSeqPartitions(ParaHypergraph &hgraph, MPI_Comm comm);
  void recursivelyBisect(const Bisection &b, MPI_Comm comm);

  void splitBisection(const Bisection &b, Bisection *&newB,
                      MPI_Comm comm) const;
  void splitBisection(const Bisection &b, Bisection *&l, Bisection *&r) const;

  int getBestPartitionProc(int cut, MPI_Comm comm) const;
  int computeMaxWt(int numBisections) const;
  double recursivelyComputeMax(double ave, int depth) const;

  inline void setNumBisectRuns(register int bRuns) { numBisectRuns = bRuns; }
};

#endif
