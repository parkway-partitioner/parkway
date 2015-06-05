
#ifndef _VCYCLEALL_BISECTION_HPP
#define _VCYCLEALL_BISECTION_HPP

// ### VCycleAllBisectionController.hpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 4/1/2005: Last Modified
//
// ###

#include "VCycleBisectionController.hpp"

using namespace std;

class VCycleAllBisectionController : public VCycleBisectionController {

protected:
public:
  VCycleAllBisectionController(int nRuns, double kT, double redFactor,
                               int eeParam, int percentile, int inc, int dispL,
                               ostream &out);
  ~VCycleAllBisectionController();

  void printVCycleType() const;
  void computeBisection();
};

#endif
