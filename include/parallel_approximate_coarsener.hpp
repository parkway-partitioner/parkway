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

#include <ostream>
#include "data_structures/bit_field.hpp"
#include "hypergraph/parallel/hypergraph.hpp"
#include "parallel_coarsener.hpp"

using parkway::hypergraph::parallel::hypergraph;
using namespace parkway::data_structures;

class parallel_approximate_coarsener : public parallel_coarsener {

protected:
  int startPercentile;
  int currPercentile;
  int increment;

public:
  parallel_approximate_coarsener(int _rank, int _numProcs, int _numParts,
                                 int percentile, int inc, std::ostream &out);

  virtual ~parallel_approximate_coarsener();
  virtual hypergraph *coarsen(hypergraph &h, MPI_Comm comm) = 0;
  virtual void set_cluster_indices(MPI_Comm comm) = 0;
  virtual void release_memory() = 0;
  virtual void display_coarsening_options() const = 0;
  virtual void build_auxiliary_structures(int numTotPins, double aveVertDeg,
                                          double aveHedgeSize) = 0;

  void load(const hypergraph &h, MPI_Comm comm);
  void compute_hyperedges_to_load(bit_field &toLoad, int numH, int *hEdgeWts,
                                  int *hEdgeOffsets, MPI_Comm comm);
};

#endif
