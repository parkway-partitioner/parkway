
#ifndef _VCYCLEFINAL_CONTROLLER_HPP
#define _VCYCLEFINAL_CONTROLLER_HPP

// ### ParaVCycleFinalController.hpp ###
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

class ParaVCycleFinalController : public ParaVCycleController {

protected:
public:
  ParaVCycleFinalController(ParaRestrCoarsener &rc, ParaCoarsener &c,
                            ParaRefiner &r, SeqController &ref, int rank,
                            int nP, int percentile, int inc, int approxRef,
                            int limit, double limitAsPercent, ostream &out);
  ~ParaVCycleFinalController();

  void runPartitioner(MPI_Comm comm);
  void printVCycleType() const;
};

#endif
