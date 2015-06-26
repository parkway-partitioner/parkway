#ifndef _PARA_REFINER_HPP
#define _PARA_REFINER_HPP
// ### ParaRefiner.hpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 4/1/2005: Last Modified
//
// ###
#include "internal/base/refiner.hpp"
#include "hypergraph/parallel/hypergraph.hpp"
#include "hypergraph/parallel/loader.hpp"
#include "data_structures/dynamic_array.hpp"
#include "data_structures/map_to_pos_int.hpp"

namespace parkway {
namespace parallel {
namespace ds = parkway::data_structures;

class refiner : public loader, public parkway::base::refiner {
 public:
  refiner(int rank, int nProcs, int nParts);

  virtual ~refiner();
  virtual void release_memory() = 0;
  virtual void initialize_data_structures(
      const parallel::hypergraph &h, MPI_Comm comm) = 0;
  virtual void reset_data_structures() = 0;
  virtual void set_partitioning_structures(int pNumber, MPI_Comm comm) = 0;
  virtual void refine(parallel::hypergraph &h, MPI_Comm comm) = 0;
  virtual int compute_cutsize(MPI_Comm comm) = 0;

  void load(const parallel::hypergraph &h, MPI_Comm comm);
  void initialize_partition_structures(const parallel::hypergraph &h,
                                       MPI_Comm comm);

  inline void set_balance_constraint(double b) {
    balance_constraint_ = b;
  }

 protected:
  int number_of_partitions_;

  ds::dynamic_array<int> partition_vector_;
  ds::dynamic_array<int> partition_vector_offsets_;
  ds::dynamic_array<int> partition_cuts_;

  int *current_partition_vector_;
  int current_partition_number_;

  double balance_constraint_;

  ds::dynamic_array<int> part_weights_;

  // newly added structures
  int number_of_non_local_vertices_;
  int number_of_non_local_vertices_to_hyperedges_;
  int *current_non_local_partition_vector_;

  ds::dynamic_array<int> non_local_vertices_;
  ds::dynamic_array<int> part_indices_;
  ds::dynamic_array<int> index_into_part_indices_;

  ds::dynamic_array<int> non_local_vertices_to_hyperedges_;
  ds::dynamic_array<int> non_local_vertices_to_hyperedges_offsets_;

  ds::map_to_pos_int to_non_local_vertices_;
};

}  // namespace parallel
}  // namespace parkway

#endif
