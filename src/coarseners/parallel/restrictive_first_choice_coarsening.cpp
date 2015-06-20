// ### ParaRestrFCCoarsener.cpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 31/12/2004: Last Modified
//
// ###

#include "coarseners/parallel/restrictive_first_choice_coarsening.hpp"

namespace parkway {
namespace parallel {

restrictive_first_choice_coarsening::restrictive_first_choice_coarsening(
    int rank, int nProcs, int nParts, int verVisOrder, int divByWt,
    int divByLen, std::ostream &out)
    : restrictive_coarsening(rank, nProcs, nParts, out) {
  vertex_visit_order_ = verVisOrder;
  divide_by_cluster_weight_ = divByWt;
  divide_by_hyperedge_length_ = divByLen;
  limit_on_index_during_corasening_ = 0;
}

restrictive_first_choice_coarsening::~restrictive_first_choice_coarsening() {}

void restrictive_first_choice_coarsening::display_options() const {
  switch (display_options_) {
  case SILENT:

    break;

  default:

    out_stream << "|--- PARA_RESTR_C:" << std::endl
               << "|- PFC:"
               << " r = " << reduction_ratio_ << " min = " << minimum_nodes_
               << " vvo = ";
      print_visit_order(vertex_visit_order_);
    out_stream << " divWt = " << divide_by_cluster_weight_ << " divLen = " <<
                                                              divide_by_hyperedge_length_
               << std::endl
               << "|" << std::endl;
    break;
  }
}

void restrictive_first_choice_coarsening::build_auxiliary_structures(int numPins,
                                                      double aveVertDeg,
                                                      double aveHedgeSize) {}

void restrictive_first_choice_coarsening::release_memory() {
  hyperedge_weights_.reserve(0);
  hyperedge_offsets_.reserve(0);
  local_pin_list_.reserve(0);

  vertex_to_hyperedges_offset_.reserve(0);
  vertex_to_hyperedges_.reserve(0);
  allocated_hyperedges_.reserve(0);

  free_memory();
}

hypergraph *restrictive_first_choice_coarsening::coarsen(hypergraph &h,
                                                         MPI_Comm comm) {
  load(h, comm);

  if (number_of_vertices_ < minimum_nodes_ || h.dont_coarsen()) {
    return nullptr;
  }

  int i;
  int j;

  int index = 0;
  int numNotMatched = number_of_local_vertices_;
  int aveVertexWt = static_cast<int>(
      ceil(static_cast<double>(total_hypergraph_weight_) / number_of_vertices_));
  int bestMatch;
  int bestMatchWt = -1;
  int endOffset1;
  int endOffset2;
  int hEdgeWt;
  int neighVertex;
  int neighVertexEntry;
  int numNeighbours;
  int vertexPart;
  int maxLocWt = 0;
  int maxWt;
  int hEdgeLen;

  int hEdge;
  int vertex;

  double metricVal;
  double maxMatchMetric;
  double reducedBy = 1;

  ds::dynamic_array<int> neighVerts;
  ds::dynamic_array<int> neighPairWts;
  ds::dynamic_array<double> connectVals;
  ds::dynamic_array<int> vertexAdjEntry(number_of_local_vertices_);
  ds::dynamic_array<int> vertices(number_of_local_vertices_);

  permute_vertices_array(vertices, number_of_local_vertices_);

  for (i = 0; i < number_of_local_vertices_; ++i)
    vertexAdjEntry[i] = -1;

  if (display_options_ > 1) {
    for (i = 0; i < number_of_local_vertices_; ++i) {
      if (vertex_weights_[i] > maxLocWt)
        maxLocWt = vertex_weights_[i];
    }

    MPI_Reduce(&maxLocWt, &maxWt, 1, MPI_INT, MPI_MAX, 0, comm);

    if (rank_ == 0) {
      out_stream << " " << maximum_vertex_weight_ << " " << maxWt << " " << aveVertexWt
                 << " ";
    }
  }

  metricVal = static_cast<double>(number_of_local_vertices_) / reduction_ratio_;

  cluster_weights_.resize(1024);
  part_vector_.resize(1024);

  limit_on_index_during_corasening_ =
      number_of_local_vertices_ - static_cast<int>(floor(metricVal - 1.0));
  cluster_index_ = 0;

  for (; index < number_of_local_vertices_; ++index) {
    if (match_vector_[vertices[index]] == -1) {
      vertex = vertices[index];
#ifdef DEBUG_COARSENER
      assert(vertex >= 0 && vertex < numLocalVertices);
#endif
      vertexPart = partition_vector_[vertex];
      endOffset1 = vertex_to_hyperedges_offset_[vertex + 1];
      numNeighbours = 0;

      for (i = vertex_to_hyperedges_offset_[vertex]; i < endOffset1; ++i) {
        hEdge = vertex_to_hyperedges_[i];
        hEdgeWt = hyperedge_weights_[hEdge];
        endOffset2 = hyperedge_offsets_[hEdge + 1];
        hEdgeLen = endOffset2 - hyperedge_offsets_[hEdge];

        for (j = hyperedge_offsets_[hEdge]; j < endOffset2; ++j) {
#ifdef DEBUG_COARSENER
          assert(locPinList[j] >= 0 && locPinList[j] < numLocalVertices);
#endif
          neighVertex = local_pin_list_[j];

          if (neighVertex != vertex &&
              partition_vector_[neighVertex] == vertexPart) {
            // ###
            // now try not to check the weight before checking the vertex
            // ###

            neighVertexEntry = vertexAdjEntry[neighVertex];

            if (neighVertexEntry == -1) {
              if (match_vector_[neighVertex] == -1)
                neighPairWts.assign(numNeighbours,
                                    vertex_weights_[vertex] + vertex_weights_[neighVertex]);
              else
                neighPairWts.assign(
                    numNeighbours,
                    vertex_weights_[vertex] +
                        cluster_weights_[match_vector_[neighVertex]]);

              neighVerts.assign(numNeighbours, neighVertex);

              if (divide_by_hyperedge_length_)
                connectVals.assign(numNeighbours, static_cast<double>(hEdgeWt) /
                                                      (hEdgeLen - 1));
              else
                connectVals.assign(numNeighbours, static_cast<double>(hEdgeWt));

              vertexAdjEntry[neighVertex] = numNeighbours;
              ++numNeighbours;
            } else {
#ifdef DEBUG_COARSENER
              assert(neighVertexEntry >= 0 && neighVertexEntry < numNeighbours);
#endif
              if (divide_by_hyperedge_length_)
                connectVals[neighVertexEntry] +=
                    (static_cast<double>(hEdgeWt) / (hEdgeLen - 1));
              else
                connectVals[neighVertexEntry] += (static_cast<double>(hEdgeWt));
            }
          }
        }
      }

      // ###
      // now choose best match from adjacent vertices
      // visited above
      // ###

      if (numNeighbours == 0) {
        // match vertex v as singleton as not connected

        match_vector_[vertex] = cluster_index_;
        part_vector_[cluster_index_] = partition_vector_[vertex];
        cluster_weights_[cluster_index_++] = vertex_weights_[vertex];
        --numNotMatched;
      } else {
        // ###
        // pick best match
        // ###

        bestMatch = -1;
        maxMatchMetric = 0.0;

        for (j = 0; j < numNeighbours; ++j) {
          if (neighPairWts[j] < maximum_vertex_weight_) {
            metricVal = connectVals[j];

            if (divide_by_cluster_weight_)
              metricVal /= neighPairWts[j];

            if (metricVal > maxMatchMetric) {
              maxMatchMetric = metricVal;
              bestMatch = neighVerts[j];
              bestMatchWt = neighPairWts[j];
            }
          }

          vertexAdjEntry[neighVerts[j]] = -1;
        }

        if (bestMatch == -1) {
          // no best match as all would violate max
          // vertex weight - match as singleton

          match_vector_[vertex] = cluster_index_;
          part_vector_[cluster_index_] = partition_vector_[vertex];
          cluster_weights_[cluster_index_++] = vertex_weights_[vertex];
          --numNotMatched;
        } else {
#ifdef DEBUG_COARSENER
          assert(bestMatch >= 0 && bestMatch < numLocalVertices);
#endif
          if (match_vector_[bestMatch] == -1) {
            // match with another unmatched

            match_vector_[vertex] = cluster_index_;
            match_vector_[bestMatch] = cluster_index_;
            part_vector_[cluster_index_] = partition_vector_[vertex];
            cluster_weights_[cluster_index_++] = bestMatchWt;
            numNotMatched -= 2;
          } else {
            // match with existing cluster of coarse vertices

            match_vector_[vertex] = match_vector_[bestMatch];
            cluster_weights_[match_vector_[vertex]] += vertex_weights_[vertex];
            --numNotMatched;
          }
        }
      }

      // check if hypergraph sufficiently shrunk

      if (index > limit_on_index_during_corasening_) {
        reducedBy = static_cast<double>(number_of_local_vertices_) /
                    (numNotMatched + cluster_index_);

        if (reducedBy > reduction_ratio_)
          break;
      }
    }
  }

  for (; index < number_of_local_vertices_; ++index) {
    vertex = vertices[index];

    if (match_vector_[vertex] == -1) {
      match_vector_[vertex] = cluster_index_;
      part_vector_[cluster_index_] = partition_vector_[vertex];
      cluster_weights_[cluster_index_++] = vertex_weights_[vertex];
    }
  }

#ifdef DEBUG_COARSENER
  for (i = 0; i < numLocalVertices; ++i) {
    assert(matchVector[i] != -1);
    assert(matchVector[i] >= 0 && matchVector[i] < clusterIndex);
  }
#endif

  set_cluster_indices(comm);

  if (static_cast<double>(number_of_vertices_) / total_clusters_ <
      MIN_ALLOWED_REDUCTION_RATIO)
    stop_coarsening_ = 1;
  else
    stop_coarsening_ = 0;

  // ###
  // now construct the coarse hypergraph using the matching vector
  // ###

  return (contract_hyperedges(h, comm));
}

void restrictive_first_choice_coarsening::permute_vertices_array(
    dynamic_array<int> &verts, int nLocVerts) {
  int i;
  switch (vertex_visit_order_) {
  case INCREASING_ORDER:
    for (i = 0; i < nLocVerts; ++i) {
      verts[i] = i;
    }
    break;

  case DECREASING_ORDER:
    for (i = 0; i < nLocVerts; ++i) {
      verts[i] = nLocVerts - i - 1;
    }
    break;

  case RANDOM_ORDER:
    for (i = 0; i < nLocVerts; ++i) {
      verts[i] = i;
    }
    verts.random_permutation();
    break;

  case INCREASING_WEIGHT_ORDER:
    for (i = 0; i < nLocVerts; ++i) {
      verts[i] = i;
    }
    verts.sort_between_using_another_array(
        0, nLocVerts - 1, vertex_weights_,
        parkway::utility::sort_order::INCREASING);
    break;

  case DECREASING_WEIGHT_ORDER:
    for (i = 0; i < nLocVerts; ++i) {
      verts[i] = i;
    }
    verts.sort_between_using_another_array(
        0, nLocVerts - 1, vertex_weights_,
        parkway::utility::sort_order::DECREASING);
    break;

  default:
    for (i = 0; i < nLocVerts; ++i) {
      verts[i] = i;
    }
    verts.random_permutation();
    break;
  }
}




void restrictive_first_choice_coarsening::set_cluster_indices(MPI_Comm comm) {
  int i;

  MPI_Scan(&cluster_index_, &minimum_cluster_index_, 1, MPI_INT, MPI_SUM, comm);

  if (rank_ == processors_ - 1)
    total_clusters_ = minimum_cluster_index_;

  minimum_cluster_index_ -= cluster_index_;

  MPI_Bcast(&total_clusters_, 1, MPI_INT, processors_ - 1, comm);

  for (i = 0; i < number_of_local_vertices_; ++i) {
#ifdef DEBUG_COARSENER
    assert(matchVector[i] != -1);
#endif
    match_vector_[i] += minimum_cluster_index_;
  }
}

void restrictive_first_choice_coarsening::print_visit_order(
    int variable) const {
  switch (variable) {
  case INCREASING_ORDER:
    out_stream << "inc-idx";
    break;

  case DECREASING_ORDER:
    out_stream << "dec-idx";
    break;

  case INCREASING_WEIGHT_ORDER:
    out_stream << "inc_wt";
    break;

  case DECREASING_WEIGHT_ORDER:
    out_stream << "dec-wt";
    break;

  default:
    out_stream << "rand";
    break;
  }
}

}  // namespace parallel
}  // namespace parkway
