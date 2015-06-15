
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

#include "parallel_controller.hpp"

using namespace std;

class basic_parallel_controller : public parallel_controller {

protected:
public:
  basic_parallel_controller(parallel_coarsener &c, parallel_refiner &r, sequential_controller &ref,
                      int rank, int nP, int percentile, int inc, int approxRef,
                      ostream &out);
  ~basic_parallel_controller();

  void reset_structures();
  void display_options() const;
  void run(MPI_Comm comm);
};

#endif
