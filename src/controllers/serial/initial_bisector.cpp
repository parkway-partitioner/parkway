// ### InitBisector.cpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 12/1/2005: Last Modified
//
// ###
#include "controllers/serial/initial_bisector.hpp"
#include "utility/logging.hpp"

namespace parkway {
namespace serial {

initial_bisector::initial_bisector(int nRuns, int insMethod, int ee)
    : fm_refiner(-1, insMethod, ee) {
  number_of_initial_runs_ = nRuns;
}

initial_bisector::~initial_bisector() {}

void initial_bisector::display_options() const {
  info("|- GIB: runs = %i\n|\n", number_of_initial_runs_);
}

void initial_bisector::initialize_bisector(serial::hypergraph &h) {
  int i;

  int gain;
  int balanced;

  h.set_number_of_partitions(number_of_initial_runs_);

  load_for_refinement(h);
  remove_ee_threshold();
  build_buckets();

  for (i = 0; i < number_of_initial_runs_; ++i) {
    set_partition_vector(i);
    partitionCutsizes[i] = set_base_vertex();
    initialise_partition_structure();
    initialise_gains_1_to_0();
    partitionCutsizes[i] -= greedy_pass();

    balanced = (part_weights_[0] <= maximum_part_weight_ && part_weights_[1] <= maximum_part_weight_);

    if (balanced)
      do {
        gain = fm_pass();
#ifdef DEBUG_REFINER
        assert(part_weights_[0] <= maximum_part_weight_);
        assert(part_weights_[1] <= maximum_part_weight_);
#endif
        partitionCutsizes[i] -= gain;
      } while (gain > 0);

    restore_buckets();
  }

  destroy_buckets();
}

int initial_bisector::set_base_vertex() {
  int i;
  int cut = 0;

  int endOffset;
  int baseVertex = RANDOM(0, numVertices);

  part_weights_[1] = 0;

  for (i = 0; i < numVertices; ++i) {
    partition_vector_[i] = 1;
    part_weights_[1] += vWeight[i];
  }

  partition_vector_[baseVertex] = 0;
  part_weights_[1] -= vWeight[baseVertex];
  part_weights_[0] = vWeight[baseVertex];

  endOffset = vOffsets[baseVertex + 1];

  for (i = vOffsets[baseVertex]; i < endOffset; ++i)
    cut += hEdgeWeight[vToHedges[i]];

  return cut;
}

int initial_bisector::choose_best_vertex_1_to_0() {
#ifdef DEBUG_REFINER
  assert((*bucketArrays[1])[maxGainEntries[1]].next);
#endif

  return ((*bucket_arrays_[1])[maximum_gain_entries_[1]].next->vertex_id);
}

int initial_bisector::greedy_pass() {
  int gain = 0;
  int bestVertex;

  while (part_weights_[1] > part_weights_[0]) {
    bestVertex = choose_best_vertex_1_to_0();

#ifdef DEBUG_REFINER
    assert(partitionVector[bestVertex] == 1);
#endif

    gain += vertex_gains_[bestVertex];

    remove_from_bucket_array(bestVertex, 1, vertex_gains_[bestVertex]);
    update_gains_1_to_0(bestVertex);
  }

  remove_buckets_from_1();

  // ###
  // if now unbalanced, do a rebalancing pass
  // ###

  if (part_weights_[0] > maximum_part_weight_)
    gain += rebalancing_pass(0);

  return gain;
}

}  // namespace serial
}  // namespace parkway
