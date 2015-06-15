#ifndef _PARA_APPROX_COARSENER_HPP
#define _PARA_APPROX_COARSENER_HPP
// ### ParaApproxCoarsener.hpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 31/12/2004: Last Modified
//
// ###
#include <iostream>
#include "data_structures/bit_field.hpp"
#include "hypergraph/parallel/hypergraph.hpp"
#include "coarseners/parallel/coarsener.hpp"

namespace parkway {
namespace parallel {
namespace ds = parkway::data_structures;

class approximate_coarsener : public coarsener {
 public:
  approximate_coarsener(int _rank, int _numProcs, int _numParts,
                        int percentile, int inc, std::ostream &out);

  virtual ~approximate_coarsener();

  virtual void set_cluster_indices(MPI_Comm comm) = 0;

  void compute_hyperedges_to_load(bit_field &toLoad, int numH, int *hEdgeWts,
                                  int *hEdgeOffsets, MPI_Comm comm);

 protected:
  int startPercentile;
  int currPercentile;
  int increment;

};

}  // namespace parallel
}  // namespace parkway

#endif
