
#ifndef _VCYCLE_CONTROLLER_HPP
#define _VCYCLE_CONTROLLER_HPP

// ### ParaVCycleController.hpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 4/1/2005: Last Modified
//
// ###

#include "ParaController.hpp"
#include "ParaRestrFCCoarsener.hpp"
#include "hypergraph/parallel/hypergraph.hpp"

using parkway::hypergraph::parallel::hypergraph;
typedef dynamic_array<int> IntArray;

class ParaVCycleController : public ParaController {

protected:
  int limitOnCycles;
  int minInterVertIndex;

  double limAsPercentOfCut;

  stack<int> numLocCurrVerts;
  stack<int> minLocCurrVertId;
  stack<IntArray *> bestVCyclePartition;

  dynamic_array<int> mapToInterVerts;

  ParaRestrCoarsener &restrCoarsener;

public:
  ParaVCycleController(ParaRestrCoarsener &rc, parallel_coarsener &c, ParaRefiner &r,
                       SeqController &ref, int rank, int nP, int percentile,
                       int inc, int approxRef, int limit, double limitAsPercent,
                       std::ostream &out);

  virtual ~ParaVCycleController();
  virtual void runPartitioner(MPI_Comm comm) = 0;
  virtual void printVCycleType() const = 0;

  void setWeightConstraints(MPI_Comm comm);
  void dispParaControllerOptions() const;

  void recordVCyclePartition(const hypergraph &h, int numIteration);
  void gatherInVCyclePartition(hypergraph &h, int cut, MPI_Comm comm);

  void projectVCyclePartition(hypergraph &cG, hypergraph &fG,
                              MPI_Comm comm);
  void shuffleVCycleVertsByPartition(hypergraph &h, MPI_Comm comm);
  // void randomVCycleVertShuffle(ParaHypergraph &h, ParaHypergraph &fineH,
  // MPI_Comm comm);
  void shiftVCycleVertsToBalance(hypergraph &h, MPI_Comm comm);
  void updateMapToOrigVerts(MPI_Comm comm);
  void resetStructs();
};

#endif
