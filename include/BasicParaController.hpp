
#ifndef _BASIC_PARA_CONTROLLER_HPP
#define _BASIC_PARA_CONTROLLER_HPP

// ### BasicParaController.hpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 4/1/2005: Last Modified
//
// ###

#include "ParaController.hpp"

using namespace std;

class BasicParaController : public ParaController {

protected:
public:
  BasicParaController(parallel_coarsener &c, ParaRefiner &r, SeqController &ref,
                      int rank, int nP, int percentile, int inc, int approxRef,
                      ostream &out);
  ~BasicParaController();

  void resetStructs();
  void dispParaControllerOptions() const;
  void runPartitioner(MPI_Comm comm);
};

#endif
