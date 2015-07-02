#ifndef _PARA_COARSENER_HPP
#define _PARA_COARSENER_HPP
// ### ParaCoarsener.hpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 4/1/2005: Last Modified
//
// ###
#include "internal/base/coarsener.hpp"

#include <limits.h>
#include "hypergraph/parallel/hypergraph.hpp"
#include "hypergraph/parallel/loader.hpp"
#include "data_structures/dynamic_array.hpp"

namespace parkway {
namespace parallel {
namespace ds = parkway::data_structures;

class coarsener : public loader, public parkway::base::coarsener {
 public:
  coarsener(int rank, int number_of_processors, int number_of_parts,
            std::ostream &out, std::string object_name = "Parallel Coarsener");

  virtual ~coarsener();

  virtual hypergraph *coarsen(hypergraph &h, MPI_Comm comm) = 0;
  virtual void set_cluster_indices(MPI_Comm comm) = 0;
  virtual void release_memory() = 0;
  virtual void display_options() const = 0;
  virtual void build_auxiliary_structures(int numTotPins, double aveVertDeg,
                                          double aveHedgeSize) = 0;

  virtual void load(const hypergraph &h, MPI_Comm comm);

  virtual hypergraph *contract_hyperedges(hypergraph &h, MPI_Comm comm);

  inline void set_balance_constraint(double constraint) {
    balance_constraint_ = constraint;
  }

  inline void set_total_hypergraph_weight(int total_weigtht) {
    total_hypergraph_weight_ = total_weigtht;
  }

 protected:
  int total_hypergraph_weight_;
  int stop_coarsening_;
  int cluster_index_;
  int total_number_of_clusters_;
  int minimum_cluster_index_;
  double balance_constraint_;

  ds::dynamic_array<int> cluster_weights_;

  void update_hypergraph_information(const hypergraph &h);
  void initialize_vertex_to_hyperedges();
  void load_non_local_hyperedges();
  void prepare_data_to_send(int n_local_hyperedges,
                            int n_local_pins,
                            int vertices_per_processor,
                            dynamic_array<int> &local_hyperedge_weights,
                            dynamic_array<int> &local_hyperedge_offsets,
                            dynamic_array<int> &local_pins,
                            MPI_Comm comm,
                            bool check_limit = false,
                            int limit = INT_MAX);

  bool within_vertex_index_range(int value) const;
};

}  // namespace parallel
}  // namespace parkway

#endif
