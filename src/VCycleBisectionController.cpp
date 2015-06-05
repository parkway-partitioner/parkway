
#ifndef _VCYCLE_BISECTION_CPP
#define _VCYCLE_BISECTION_CPP

// ### VCycleBisectionController.cpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 4/1/2005: Last Modified
//
// ###

#include "VCycleBisectionController.hpp"

VCycleBisectionController::VCycleBisectionController(
    const int nRuns, const double kT, const double redFactor, int eeParam,
    int percentile, int inc, int dispL, ostream &out)
    : BisectionController(nRuns, kT, redFactor, eeParam, percentile, inc, dispL,
                          out) {
  restrCoarsener = NULL;

  vCyclePartition.setLength(0);
}

VCycleBisectionController::~VCycleBisectionController() {
  DynaMem<RestrCoarsener>::deletePtr(restrCoarsener);
}

void VCycleBisectionController::dispBisectionControllerOptions() const {
  switch (dispLevel) {
  case SILENT:
    break;

  default:

    out_stream << "|- V-CYC BSECTOR:";
#ifdef DEBUG_CONTROLLER
    assert(coarsener && initBisector && refiner && restrCoarsener);
#endif
    out_stream << " kT = " << keepThreshold << " rF = " << reductionFactor
               << " %le = " << startPercentile
               << " %inc = " << percentileIncrement;
    printVCycleType();
    out_stream << endl << "|" << endl;
    coarsener->dispCoarsenerOptions(out_stream);
    restrCoarsener->dispCoarsenerOptions(out_stream);
    initBisector->dispInitBisectorOptions(out_stream);
    refiner->dispRefinerOptions(out_stream);

    break;
  }
}

void VCycleBisectionController::buildRestrCoarsener(double redRatio, int cType,
                                                    int minNodes) {
  switch (cType) {
  case RestrFCwithFanOutDiv:
    restrCoarsener = new RestrFCCoarsener(Shiftl(minNodes, 1), -1, redRatio, 1,
                                          1, dispLevel);
    break;

  case RestrFCwithoutFanOutDiv:
    restrCoarsener = new RestrFCCoarsener(Shiftl(minNodes, 1), -1, redRatio, 0,
                                          1, dispLevel);
    break;

  case RestrFCwithFanOut:
    restrCoarsener = new RestrFCCoarsener(Shiftl(minNodes, 1), -1, redRatio, 1,
                                          0, dispLevel);
    break;

  case RestrFCwithoutFanOut:
    restrCoarsener = new RestrFCCoarsener(Shiftl(minNodes, 1), -1, redRatio, 0,
                                          0, dispLevel);
    break;

  default:
    restrCoarsener = new RestrFCCoarsener(Shiftl(minNodes, 1), -1, redRatio, 1,
                                          1, dispLevel);
    break;
  }
}

void VCycleBisectionController::recordVCyclePartition(
    const int *pVector, int numV) {
  int i;

  for (i = 0; i < numV; ++i) {
#ifdef DEBUG_CONTROLLER
    assert(pVector[i] >= 0 && pVector[i] < 2);
#endif
    vCyclePartition[i] = pVector[i];
  }
}

void VCycleBisectionController::storeBestPartition(const int *pVector,
                                                   int numV) {
  int i;

  for (i = 0; i < numV; ++i) {
#ifdef DEBUG_CONTROLLER
    assert(pVector[i] >= 0 && pVector[i] < 2);
#endif
    bestPartition[i] = pVector[i];
  }
}

#endif
