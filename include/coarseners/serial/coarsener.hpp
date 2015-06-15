#ifndef _COARSENER_HPP
#define _COARSENER_HPP
// ### Coarsener.hpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 4/1/2005: Last Modified
//
// ###
#include "coarseners/base_coarsener.hpp"
#include "hypergraph/serial/hypergraph.hpp"
#include "hypergraph/serial/loader.hpp"

namespace parkway {
namespace serial {

class coarsener : public loader, public parkway::coarsener::base_coarsener {
 public:
  coarsener(int min, int maxwt, double ratio, int dispL,
            std::string object_name = "Serial Coarsener");
  virtual ~coarsener();
  virtual hypergraph *coarsen(const hypergraph &h) = 0;
  virtual void display_options(std::ostream &out) const = 0;

  hypergraph *build_coarse_hypergraph(int *coarse_weights,
                                      int number_of_coarse_vertices,
                                      int total_weight) const;
};

}  // namespace serial
}  // namespace parkway

#endif
