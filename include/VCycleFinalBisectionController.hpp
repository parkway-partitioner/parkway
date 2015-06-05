
#ifndef _VCYCLEFINAL_BISECTION_HPP
#define _VCYCLEFINAL_BISECTION_HPP

// ### VCycleFinalBisectionController.hpp ###
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

class VCycleFinalBisectionController : public VCycleBisectionController {

protected:
public:
  VCycleFinalBisectionController(int nRuns, double kT, double redFactor,
                                 int eeParam, int percentile, int inc,
                                 int dispL, ostream &out);
  ~VCycleFinalBisectionController();

  void printVCycleType() const;
  void computeBisection();
};

#endif
