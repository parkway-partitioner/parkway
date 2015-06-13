
#ifndef _VCYCLEALL_CONTROLLER_HPP
#define _VCYCLEALL_CONTROLLER_HPP

// ### ParaVCycleController.hpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 4/1/2005: Last Modified
//
// ###

#include "ParaVCycleController.hpp"

using namespace std;

class ParaVCycleAllController : public ParaVCycleController {

protected:
public:
  ParaVCycleAllController(ParaRestrCoarsener &rc, parallel_coarsener &c,
                          ParaRefiner &r, SeqController &ref, int rank, int nP,
                          int percentile, int inc, int approxRef, int limit,
                          double limitAsPercent, ostream &out);
  ~ParaVCycleAllController();

  void runPartitioner(MPI_Comm comm);
  void printVCycleType() const;
};

#endif
