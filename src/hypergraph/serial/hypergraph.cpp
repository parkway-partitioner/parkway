#ifndef _HYPERGRAPH_CPP
#define _HYPERGRAPH_CPP

// ### Hypergraph.cpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 5/1/2005: Last Modified
//
// ###

#include "hypergraph/serial/hypergraph.hpp"
#include <cmath>

namespace parkway {
namespace serial {

hypergraph::hypergraph(int *vWts, int numVerts)
    : hg::base_hypergraph(numVerts) {

  match_vector_.reserve(number_of_vertices_);
  vertex_weights_.set_data(vWts, number_of_vertices_);

  for (int i = 0; i < number_of_vertices_; ++i) {
    match_vector_[i] = -1;
  }
}

hypergraph::hypergraph(int *vWts, int *p_vector, int numVerts, int cut)
    : hypergraph(vWts, numVerts) {

  partition_vector_.set_data(p_vector, number_of_vertices_);

  partition_cuts_.reserve(1);
  partition_vector_offsets_.reserve(2);

  partition_vector_offsets_[0] = 0;
  partition_vector_offsets_[1] = number_of_vertices_;
  partition_cuts_[0] = cut;
}

hypergraph::~hypergraph() {}

void hypergraph::buildVtoHedges() {
#ifdef DEBUG_HYPERGRAPH
  assert(numVertices >= 0);
  assert(numHedges >= 0);
  assert(numPins >= 0);
  assert(pinList.getLength() >= 0);
  assert(hEdgeOffsets.getLength() >= 0);
#endif

  ds::dynamic_array<int> vDegs(number_of_vertices_);
  vertex_to_hyperedges_.reserve(number_of_pins_);
  vertex_offsets_.reserve(number_of_vertices_ + 1);

  for (int i = 0; i < number_of_vertices_; ++i) {
    vDegs[i] = 0;
  }

  for (int i = 0; i < number_of_pins_; ++i) {
#ifdef DEBUG_HYPERGRAPH
    assert(pinList[i] >= 0 && pinList[i] < numVertices);
#endif
    ++vDegs[pin_list_[i]];
  }

  int i;
  int j = 0;
  for (i = 0; i < number_of_vertices_; ++i) {
    vertex_offsets_[i] = j;
    j += vDegs[i];
    vDegs[i] = 0;
  }
  vertex_offsets_[i] = j;

  for (i = 0; i < number_of_hyperedges_; ++i) {
    int endOffset = hyperedge_offsets_[i + 1];
    for (j = hyperedge_offsets_[i]; j < endOffset; ++j) {
      int ij = pin_list_[j];
      vertex_to_hyperedges_[vertex_offsets_[ij] + (vDegs[ij]++)] = i;
    }
  }
}

void hypergraph::reset_match_vector() {
  for (int i = 0; i < number_of_vertices_; ++i) {
    match_vector_[i] = -1;
  }
}

void hypergraph::reset_partition_vector() {
  partition_vector_.reserve(0);
  partition_cuts_.reserve(0);
  partition_vector_offsets_.reserve(0);
  number_of_partitions_ = 0;
}

void hypergraph::reset_vertex_maps() {
  reset_match_vector();
  reset_partition_vector();
}

void hypergraph::project_partitions(const hypergraph &coarse_graph) {
  number_of_partitions_ = coarse_graph.number_of_partitions();

  int *coarse_partition_vector = coarse_graph.partition_vector();
  int *coarse_partition_offsets = coarse_graph.partition_offsets();
  int *coarse_partition_cuts = coarse_graph.partition_cuts();

  partition_cuts_.reserve(number_of_partitions_);
  partition_vector_offsets_.reserve(number_of_partitions_ + 1);
  partition_vector_.reserve(number_of_partitions_ * number_of_vertices_);

  int j = 0;
  for (int i = 0; i <= number_of_partitions_; ++i) {
    partition_vector_offsets_[i] = j;
    j += number_of_vertices_;
  }

  for (int i = 0; i < number_of_partitions_; ++i) {
    partition_cuts_[i] = coarse_partition_cuts[i];

    int index = coarse_partition_offsets[i];
    int start_offset = partition_vector_offsets_[i];

    for (j = 0; j < number_of_vertices_; ++j) {
      int v = match_vector_[j];
      partition_vector_[start_offset + j] = coarse_partition_vector[index + v];
    }
  }
}

void hypergraph::remove_bad_partitions(double fraction_ok) {
  int best_cut = partition_cuts_[0];
  int new_partitions = 0;

  // Find the best cut.
  for (int i = 1; i < number_of_partitions_; ++i) {
    if (partition_cuts_[i] < best_cut) {
      best_cut = partition_cuts_[i];
    }
  }

  int accepted_cut = best_cut + static_cast<int>(
      std::floor(static_cast<double>(best_cut) * fraction_ok));

  int index_into_new = 0;
  int index_into_old = 0;

  for (int i = 0; i < number_of_partitions_; ++i) {
    if (partition_cuts_[i] <= accepted_cut) {
      bool seen_before = false;

      for (int j = 0; j < new_partitions; ++j) {
        if (partition_cuts_[j] == partition_cuts_[i]) {
          seen_before = true;
        }
      }

      if (!seen_before) {
        if (index_into_old > index_into_new) {
          int end_offset = partition_vector_offsets_[i + 1];

          while (index_into_old < end_offset) {
            partition_vector_[index_into_new++] =
                  partition_vector_[index_into_old];
            ++index_into_old;
          }
        } else {
          index_into_old += number_of_vertices_;
          index_into_new += number_of_vertices_;
        }

        partition_cuts_[new_partitions++] = partition_cuts_[i];
      } else {
        index_into_old += number_of_vertices_;
      }
    } else {
      index_into_old += number_of_vertices_;
    }
  }

  number_of_partitions_ = new_partitions;
}

void hypergraph::set_number_of_partitions(int nPartitions) {
  number_of_partitions_ = nPartitions;

  partition_cuts_.reserve(number_of_partitions_);
  partition_vector_.reserve(number_of_partitions_ * number_of_vertices_);
  partition_vector_offsets_.reserve(number_of_partitions_ + 1);

  int j = 0;
  for (int i = 0; i <= number_of_partitions_; ++i) {
    partition_vector_offsets_[i] = j;
    j += number_of_vertices_;
  }
}

void hypergraph::copy_out_partition(int *p_vector, int nV, int pNo) const {
#ifdef DEBUG_HYPERGRAPH
  assert(pNo >= 0 && pNo < numPartitions);
  assert(nV == numVertices);
#endif

  int *vOffset = &partition_vector_[partition_vector_offsets_[pNo]];
  for (int i = 0; i < number_of_vertices_; ++i) {
    p_vector[i] = vOffset[i];
  }
}

void hypergraph::copy_in_partition(const int *p_vector, int nV, int pNo,
                                   int cut) {
#ifdef DEBUG_HYPERGRAPH
  assert(pNo >= 0 && pNo < numPartitions);
  assert(nV == numVertices);
#endif

  int *vOffset = &partition_vector_[partition_vector_offsets_[pNo]];
  for (int i = 0; i < number_of_vertices_; ++i) {
    vOffset[i] = p_vector[i];
  }

  partition_cuts_[pNo] = cut;
}

void hypergraph::print_characteristics(std::ostream &out_stream) {
  out_stream << " |cGraph| " << number_of_vertices_ << " "
      << number_of_hyperedges_ << " " << number_of_pins_ << " : ";

  double weighted_ave = 0;
  double percentile_75;
  double percentile_50;
  double percentile_25;
  double percentile_95;

  ds::dynamic_array<int> hEdgeLens(number_of_hyperedges_);
  ds::dynamic_array<int> hEdges(number_of_hyperedges_);
  ds::dynamic_array<int> vertices(number_of_vertices_);

  int j = 0;
  for (int i = 0; i < number_of_hyperedges_; ++i) {
    hEdges[i] = i;
    hEdgeLens[i] = hyperedge_offsets_[i + 1] - hyperedge_offsets_[i];
    weighted_ave += (hEdgeLens[i] * hyperedge_weights_[i]);
    j += hyperedge_weights_[i];
  }

  percentile_75 = (static_cast<double>(j) * 75) / 100;
  percentile_50 = (static_cast<double>(j) * 50) / 100;
  percentile_95 = (static_cast<double>(j) * 95) / 100;
  percentile_25 = (static_cast<double>(j) * 25) / 100;

  out_stream << weighted_ave / j << " ";

  Funct::qsortByAnotherArray(0, number_of_hyperedges_ - 1, hEdges.data(),
                             hEdgeLens.data(), INC);

  int ij = 0;
  for (int i = 0; i < number_of_hyperedges_;) {
    j += hyperedge_weights_[hEdges[i++]];

    if (ij == 0 && j > percentile_25) {
      out_stream << hEdgeLens[hEdges[i]] << " ";
      ++ij;
    }

    if (ij == 1 && j > percentile_50) {
      out_stream << hEdgeLens[hEdges[i]] << " ";
      ++ij;
    }

    if (ij == 2 && j > percentile_75) {
      out_stream << hEdgeLens[hEdges[i]] << " ";
      ++ij;
    }

    if (ij == 3 && j > percentile_95) {
      out_stream << hEdgeLens[hEdges[i]] << " ";
      ++ij;
    }

    if (i == number_of_hyperedges_ - 1) {
      out_stream << hEdgeLens[hEdges[i]] << " : ";
    }
  }

  j = 0;
  for (int i = 0; i < number_of_vertices_; ++i) {
    vertices[i] = i;
    j += vertex_weights_[i];
  }

  percentile_75 = (static_cast<double>(j) * 75) / 100;
  percentile_50 = (static_cast<double>(j) * 50) / 100;
  percentile_95 = (static_cast<double>(j) * 95) / 100;
  percentile_25 = (static_cast<double>(j) * 25) / 100;

  Funct::qsortByAnotherArray(0, number_of_vertices_ - 1, vertices.data(),
                             vertex_weights_.data(), INC);

  j = 0;
  ij = 0;
  for (int i = 0; i < number_of_vertices_;) {
    j += vertex_weights_[vertices[i++]];

    if (ij == 0 && j > percentile_25) {
      out_stream << vertex_weights_[vertices[i]] << " ";
      ++ij;
    }

    if (ij == 1 && j > percentile_50) {
      out_stream << vertex_weights_[vertices[i]] << " ";
      ++ij;
    }

    if (ij == 2 && j > percentile_75) {
      out_stream << vertex_weights_[vertices[i]] << " ";
      ++ij;
    }

    if (ij == 3 && j > percentile_95) {
      out_stream << vertex_weights_[vertices[i]] << " ";
      ++ij;
    }

    if (i == number_of_vertices_ - 1) {
      out_stream << vertex_weights_[vertices[i]] << std::endl;
    }
  }
}

int hypergraph::keep_best_partition() {
  int bestPartition = 0;
  int best_cut = partition_cuts_[0];

  for (int i = 1; i < number_of_partitions_; ++i) {
    if (partition_cuts_[i] < best_cut) {
      bestPartition = i;
      best_cut = partition_cuts_[i];
    }
  }

  if (bestPartition != 0) {
    int bestOffset = partition_vector_offsets_[bestPartition];
    for (int i = 0; i < number_of_vertices_; ++i) {
      partition_vector_[i] = partition_vector_[bestOffset + i];
    }
  }

  partition_cuts_[0] = best_cut;
  number_of_partitions_ = 1;

  return best_cut;
}

int hypergraph::export_hyperedge_weight() const {
  int ij = 0;
  for (int i = 0; i < number_of_hyperedges_; ++i) {
    ij += hyperedge_weights_[i];
  }

  return ij;
}

int hypergraph::cut_size(int nP, int partitionNo) const {
  int *p_vector = &partition_vector_[partition_vector_offsets_[partitionNo]];
  ds::dynamic_array<int> spanned(nP);
  int k_1_cut = 0;

  for (int i = 0; i < number_of_hyperedges_; ++i) {
    int endOffset = hyperedge_offsets_[i + 1];
    int number_spanned = 0;

    for (int j = 0; j < nP; ++j) {
      spanned[j] = 0;
    }

    for (int j = hyperedge_offsets_[i]; j < endOffset; ++j) {
      int vertex_part = p_vector[pin_list_[j]];
#ifdef DEBUG_HYPERGRAPH
      assert(vertex_part >= 0 && vertex_part < nP);
#endif
      if (!spanned[vertex_part]) {
        spanned[vertex_part] = 1;
        ++number_spanned;
      }
    }

    k_1_cut += ((number_spanned - 1) * hyperedge_weights_[i]);
  }

  return k_1_cut;
}

int hypergraph::sum_of_external_degrees(int nP, int partitionNo) const {
  int *p_vector = &partition_vector_[partition_vector_offsets_[partitionNo]];

  ds::dynamic_array<int> spanned(nP);

  int soed = 0;

  for (int i = 0; i < number_of_hyperedges_; ++i) {
    int endOffset = hyperedge_offsets_[i + 1];
    int number_spanned = 0;

    for (int j = 0; j < nP; ++j) {
      spanned[j] = 0;
    }

    for (int j = hyperedge_offsets_[i]; j < endOffset; ++j) {
      int vertex_part = p_vector[pin_list_[j]];
#ifdef DEBUG_HYPERGRAPH
      assert(vertex_part >= 0 && vertex_part < nP);
#endif
      if (!spanned[vertex_part]) {
        spanned[vertex_part] = 1;
        ++number_spanned;
      }
    }

    if (number_spanned > 1) {
      soed += (number_spanned * hyperedge_weights_[i]);
    }
  }

  return soed;
}

void hypergraph::initialize_cut_sizes(int numParts) {
#ifdef DEBUG_HYPERGRAPH
  assert(numPartitions >= 1);
  assert(partitionCuts.getLength() >= 1);
  assert(partitionVector.getLength() >= numVertices);
  assert(partitionVectorOffsets.getLength() > 1);
#endif

  for (int i = 0; i < number_of_partitions_; ++i) {
    partition_cuts_[i] = cut_size(numParts, i);
  }
}

void hypergraph::check_partitions(int nP, int maxWt) const {
  ds::dynamic_array<int> partWts(nP);
  for (int i = 0; i < number_of_partitions_; ++i) {
    int cut = cut_size(nP, i);
    assert(partition_cuts_[i] == cut);

    for (int j = 0; j < nP; ++j) {
      partWts[j] = 0;
    }

    int *p_vector = &partition_vector_[partition_vector_offsets_[i]];

    for (int j = 0; j < number_of_vertices_; ++j) {
      partWts[p_vector[j]] += vertex_weights_[j];
    }

    check_part_weights_are_less_than(partWts.data(), nP, maxWt);
  }
}

void hypergraph::check_partition(int numPartition, int nP, int maxWt) const {
  ds::dynamic_array<int> partWts(nP);

  int cut = cut_size(nP, numPartition);
  assert(cut == partition_cuts_[numPartition]);

  for (int i = 0; i < nP; ++i)
    partWts[i] = 0;

  int *p_vector = &partition_vector_[partition_vector_offsets_[numPartition]];

  for (int i = 0; i < number_of_vertices_; ++i)
    partWts[p_vector[i]] += vertex_weights_[i];

  check_part_weights_are_less_than(partWts.data(), nP, maxWt);
}

void hypergraph::check_part_weights_are_less_than(int *part_weights,
                                                  const int number,
                                                  int maximum) const {
  for (int i = 0; i < number; ++i) {
    assert(part_weights[i] <= maximum);
  }
}

void hypergraph::convert_to_DOMACS_graph_file(const char *fN) const {
  int i;
  int j;
  int ij;
  int endOffset;
  int start_offset;
  int v1;
  int v2;
  int numEdges = 0;
  int maxVdegree = 0;

  char writeString[2064];

  std::ofstream out_stream;

  ds::dynamic_array<int> numVNeighs(number_of_vertices_);
  ds::dynamic_array<int> vDegs(number_of_vertices_);

  ds::dynamic_array<ds::dynamic_array<int> *> vNeighs(number_of_vertices_);

  for (i = 0; i < number_of_vertices_; ++i) {
    numVNeighs[i] = 0;
    vNeighs[i] = new ds::dynamic_array<int>(32);
  }

  // ###
  // first qsort the hyperedges
  // then create the adjacency list
  // ###

  for (i = 0; i < number_of_hyperedges_; ++i) {
    start_offset = hyperedge_offsets_[i];
    endOffset = hyperedge_offsets_[i + 1] - 1;

    Funct::qsort(start_offset, endOffset, pin_list_.data());

    for (j = start_offset; j < endOffset; ++j) {
      v1 = pin_list_[j];

      for (ij = j + 1; ij < endOffset; ++ij) {
        v2 = pin_list_[ij];

        if (Funct::search(vNeighs[v1]->data(), numVNeighs[v1], v2) == -1) {
          vNeighs[v1]->assign(numVNeighs[v1]++, v2);
          ++numEdges;
        }
      }
    }
  }

  for (i = 0; i < number_of_vertices_; ++i)
    vDegs[i] = 0;

  for (i = 0; i < number_of_vertices_; ++i) {
    for (j = 0; j < numVNeighs[i]; ++j) {
      v1 = (*vNeighs[i])[j];
      ++vDegs[i];
      ++vDegs[v1];
    }
  }

  for (i = 0; i < number_of_vertices_; ++i)
    if (vDegs[i] > maxVdegree)
      maxVdegree = vDegs[i];

  // ###
  // now translate the adjacency list to file
  // ###

  out_stream.open(fN, std::ofstream::out | std::ofstream::binary);

  if (!out_stream.is_open()) {
    std::cout << "error opening " << fN << std::endl;
    return;
  }

  strcpy(writeString, "c Hypergraph to Graph representation\n");
  out_stream.write(&writeString[0], strlen(writeString));

  sprintf(writeString, "p edge %d %d\n", number_of_vertices_, numEdges);
  out_stream.write(&writeString[0], strlen(writeString));

  for (i = 0; i < number_of_vertices_; ++i) {
    endOffset = numVNeighs[i];

    for (j = 0; j < endOffset; ++j) {
      v2 = (*vNeighs[i])[j];

      sprintf(writeString, "e %d %d\n", i + 1, v2 + 1);
      out_stream.write(&writeString[0], strlen(writeString));
    }
  }

  out_stream.close();

  for (i = 0; i < number_of_vertices_; ++i)
    dynamic_memory::delete_pointer<ds::dynamic_array<int> >(vNeighs[i]);
}

void hypergraph::print_percentiles(std::ostream &out_stream) {
  int i;
  int j;
  int ij;

  out_stream << " #vertices: " << number_of_vertices_
      << " #hyperedges: " << number_of_hyperedges_
      << " #pins: " << number_of_pins_ << std::endl;

  double weighted_ave = 0;
  double percentile_95;
  double percentile_75;
  double percentile_50;
  double percentile_25;

  ds::dynamic_array<int> indices;
  ds::dynamic_array<int> hEdgeLens(number_of_hyperedges_);

  /* display hyperedge information */

  j = 0;
  indices.reserve(number_of_hyperedges_);

  for (i = 0; i < number_of_hyperedges_; ++i) {
    indices[i] = i;
    j += hyperedge_weights_[i];
    hEdgeLens[i] = hyperedge_offsets_[i + 1] - hyperedge_offsets_[i];
    weighted_ave += (hEdgeLens[i] * hyperedge_weights_[i]);
  }

  percentile_95 = (static_cast<double>(j) * 95) / 100;
  percentile_75 = (static_cast<double>(j) * 75) / 100;
  percentile_50 = (static_cast<double>(j) * 50) / 100;
  percentile_25 = (static_cast<double>(j) * 25) / 100;

  out_stream << "hyperedge capacity_ percentiles: (weighted ave, 25, 50, 75, "
      "95, maxLength) " << std::endl;
  out_stream << "\t" << weighted_ave / j << " ";

  Funct::qsortByAnotherArray(0, number_of_hyperedges_ - 1, indices.data(),
                             hEdgeLens.data(), INC);

  j = 0;
  i = 0;
  ij = 0;

  for (; i < number_of_hyperedges_;) {
    j += hyperedge_weights_[indices[i++]];

    if (ij == 0 && j > percentile_25) {
      out_stream << hEdgeLens[indices[i]] << " ";
      ++ij;
    }

    if (ij == 1 && j > percentile_50) {
      out_stream << hEdgeLens[indices[i]] << " ";
      ++ij;
    }

    if (ij == 2 && j > percentile_75) {
      out_stream << hEdgeLens[indices[i]] << " ";
      ++ij;
    }

    if (ij == 3 && j > percentile_95) {
      out_stream << hEdgeLens[indices[i]] << " ";
      ++ij;
    }

    if (i == number_of_hyperedges_ - 1) {
      out_stream << hEdgeLens[indices[i]] << std::endl;
    }
  }

  /* display vertex information */

  j = 0;
  indices.reserve(number_of_vertices_);

  for (i = 0; i < number_of_vertices_; ++i) {
    indices[i] = i;
    j += vertex_weights_[i];
  }

  percentile_95 = (static_cast<double>(j) * 95) / 100;
  percentile_75 = (static_cast<double>(j) * 75) / 100;
  percentile_50 = (static_cast<double>(j) * 50) / 100;
  percentile_25 = (static_cast<double>(j) * 25) / 100;

  Funct::qsortByAnotherArray(0, number_of_vertices_ - 1, indices.data(),
                             vertex_weights_.data(), INC);

  out_stream << "vertex weight percentiles: (ave, 25, 50, 75, 95, maxWeight) "
      << std::endl;
  out_stream << "\t" << static_cast<double>(j) / number_of_vertices_ << " ";

  j = 0;
  i = 0;
  ij = 0;

  for (; i < number_of_vertices_;) {
    j += vertex_weights_[indices[i++]];

    if (ij == 0 && j > percentile_25) {
      out_stream << vertex_weights_[indices[i]] << " ";
      ++ij;
    }

    if (ij == 1 && j > percentile_50) {
      out_stream << vertex_weights_[indices[i]] << " ";
      ++ij;
    }

    if (ij == 2 && j > percentile_75) {
      out_stream << vertex_weights_[indices[i]] << " ";
      ++ij;
    }

    if (ij == 3 && j > percentile_95) {
      out_stream << vertex_weights_[indices[i]] << " ";
      ++ij;
    }

    if (i == number_of_vertices_ - 1) {
      out_stream << vertex_weights_[indices[i]] << std::endl;
    }
  }
}

}  // serial
}  // parkway

#endif
