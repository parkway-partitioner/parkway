
#ifndef _VCYCLE_BISECTION_HPP
#define _VCYCLE_BISECTION_HPP

// ### VCycleBisectionController.hpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 4/1/2005: Last Modified
//
// ###

#include "BisectionController.hpp"
#include "RestrFCCoarsener.hpp"

using namespace std;

class VCycleBisectionController : public BisectionController {
protected:
  FastDynaArray<int> vCyclePartition;

  RestrCoarsener *restrCoarsener;

public:
  VCycleBisectionController(const int nRuns, const double kT,
                            const double redFactor, int eeParam, int percentile,
                            int inc, int dispL, ostream &out);
  virtual ~VCycleBisectionController();

  void dispBisectionControllerOptions() const;
  void buildRestrCoarsener(double rRatio, int cType, int minV);
  void recordVCyclePartition(const int *pVector, int numV);
  void storeBestPartition(const int *pVector, int numV);

  virtual void printVCycleType() const = 0;
  virtual void computeBisection() = 0;
};

#endif
