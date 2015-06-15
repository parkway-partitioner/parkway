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

#include "internal/parallel_controller.hpp"
#include "internal/serial_controller.hpp"

namespace parkway {
namespace parallel {

class basic_contoller : public controller {
protected:
 public:
  basic_contoller(parallel_coarsener &c, refiner &r,
                  serial::controller &ref, int rank, int nP, int percentile,
                  int inc, int approxRef, std::ostream &out);
  ~basic_contoller();

  void reset_structures();
  void display_options() const;
  void run(MPI_Comm comm);
};

}  // namespace parallel
}  // namespace parkway

#endif
