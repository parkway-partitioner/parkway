
#ifndef _BISECTION_CONTROLLER_HPP
#define _BISECTION_CONTROLLER_HPP

// ### BisectionController.hpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 4/1/2005: Last Modified
//
// ###

#include <ostream>
#include "data_structures/stack.hpp"
#include "FCCoarsener.hpp"
#include "InitBisector.hpp"

using namespace parkway::data_structures;

class BisectionController {

protected:
  std::ostream &out_stream;

  int numSeqRuns;
  int eeParam;
  int dispLevel;

  int startPercentile;
  int percentileIncrement;

  int numOrigVertices;
  int maxPartWt;

  double keepThreshold;
  double reductionFactor;

  Coarsener *coarsener;
  Refiner *refiner;
  InitBisector *initBisector;
  stack<serial_hypergraph *> hGraphs;

  dynamic_array<int> bestPartition;

public:
  BisectionController(int nRuns, double kT, double redFactor, int eeParam,
                      int percentile, int inc, int dispL, std::ostream &out);

  virtual ~BisectionController();
  virtual void dispBisectionControllerOptions() const;

  void buildCoarsener(double redRatio, int cType, int minNodes);
  void buildInitBisector(int numInitRuns);
  void buildRefiner(int queueD);

  virtual void computeBisection();

  void bisect(serial_hypergraph *h, int maxPartWt);

  inline void setNumRuns(int r) { numSeqRuns = r; }
};

#endif
