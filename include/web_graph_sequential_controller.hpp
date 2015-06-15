

#ifndef _WEBGRAPH_SEQ_CONTROLLER_HPP
#define _WEBGRAPH_SEQ_CONTROLLER_HPP

// ### WebGraphSeqController.hpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 07/04/2005: Last Modified
//
// ###

#include "sequential_controller.hpp"
#include "hypergraph/parallel/hypergraph.hpp"

namespace parallel = parkway::parallel;

class web_graph_sequential_controller : public sequential_controller {
protected:
public:
  web_graph_sequential_controller(int rank, int nProcs, int nParts, ostream &out);

  virtual ~web_graph_sequential_controller();
  virtual void display_options() const = 0;
  virtual void run(parallel::hypergraph &hgraph, MPI_Comm comm) = 0;
  virtual void initialize_sequential_partitions(parallel::hypergraph &h,
                                                MPI_Comm comm);

  void initialize_coarsest_hypergraph(parallel::hypergraph &hgraph,
                                      MPI_Comm comm);
};

#endif
