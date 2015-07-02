// ### GreedyKwayRefiner.cpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 30/11/2004: Last Modified
//
// ###
#include "refiners/serial/greedy_k_way_refiner.hpp"
#include "utility/logging.hpp"

namespace parkway {
namespace serial {

greedy_k_way_refiner::greedy_k_way_refiner(int max, int nparts, double ave,
                                           double lim) {
  maximum_part_weight_ = max;
  number_of_parts_ = nparts;
  average_part_weight_ = ave;
  limit_ = lim;
  number_of_non_positive_moves_ = 0;
  part_weights_.resize(number_of_parts_);

  number_of_neighboring_parts_.resize(0);
  neighbors_of_vertex_.resize(0);
  neighbors_of_vertex_offsets_.resize(0);
  hyperedge_vertices_in_part_.resize(0);
  hyperedge_vertices_in_part_offsets_.resize(0);
  vertices_.resize(0);
  vertex_seen_.resize(0);
  seen_vertices_.resize(0);
  parts_spanned_.resize(0);

#ifdef DEBUG_REFINER
  assert(limit >= 0 && limit <= 1.0);
#endif
}

greedy_k_way_refiner::~greedy_k_way_refiner() {}

void greedy_k_way_refiner::display_options() const {
  info("|- GKWAY: lim = %.2f\n|\n", limit_);
}

void greedy_k_way_refiner::build_data_structures() {
  int i;

  parts_spanned_.resize(number_of_parts_);
  vertices_.resize(numVertices);
  vertex_seen_.resize(numVertices);
  seen_vertices_.resize(numVertices);
  number_of_neighboring_parts_.resize(numVertices);
  neighbors_of_vertex_offsets_.resize(numVertices + 1);
  neighbors_of_vertex_.resize(numVertices * number_of_parts_);

  hyperedge_vertices_in_part_offsets_.resize(numHedges + 1);
  hyperedge_vertices_in_part_.resize(numHedges * number_of_parts_);

  neighbors_of_vertex_offsets_[0] = 0;
  hyperedge_vertices_in_part_offsets_[0] = 0;

  for (i = 1; i <= numVertices; ++i)
    neighbors_of_vertex_offsets_[i] = neighbors_of_vertex_offsets_[i - 1] +
                                      number_of_parts_;

  for (i = 1; i <= numHedges; ++i)
    hyperedge_vertices_in_part_offsets_[i] = hyperedge_vertices_in_part_offsets_[i - 1] +
                                             number_of_parts_;
}

void greedy_k_way_refiner::destroy_data_structures() {
  parts_spanned_.resize(0);
  vertices_.resize(0);
  vertex_seen_.resize(0);
  seen_vertices_.resize(0);
  number_of_neighboring_parts_.resize(0);
  neighbors_of_vertex_offsets_.resize(0);
  neighbors_of_vertex_.resize(0);

  hyperedge_vertices_in_part_offsets_.resize(0);
  hyperedge_vertices_in_part_.resize(0);
}

int greedy_k_way_refiner::initialize_data_structures() {
  int i;
  int j;
  int ij;

  int v;
  int endIndex;
  int vertOffset;
  int partOffset;
  int numPartsSpanned;
  int vPart;
  int k_1cut = 0;

  for (i = 0; i < number_of_parts_; ++i)
    part_weights_[i] = 0;

  for (i = 0; i < numVertices; ++i) {
#ifdef DEBUG_REFINER
    assert(partitionVector[i] >= 0 && partitionVector[i] < numParts);
#endif
    part_weights_[partition_vector_[i]] += vWeight[i];
    number_of_neighboring_parts_[i] = 0;
    vertex_seen_[i] = -1;
  }

  endIndex = neighbors_of_vertex_offsets_[numVertices];

  for (i = 0; i < endIndex; ++i)
    neighbors_of_vertex_[i] = 0;

  endIndex = hyperedge_vertices_in_part_offsets_[numHedges];

  for (i = 0; i < endIndex; ++i)
    hyperedge_vertices_in_part_[i] = 0;

  // ###
  // initialise the hyperedge structures
  // and the vertices structure
  // ###

  for (i = 0; i < numHedges; ++i) {
    endIndex = hEdgeOffsets[i + 1];
    numPartsSpanned = 0;

    // ###
    // update the hyperedge part distributions
    // ###

    for (j = hEdgeOffsets[i]; j < endIndex; ++j) {
      vPart = partition_vector_[pinList[j]];
      partOffset = hyperedge_vertices_in_part_offsets_[i] + vPart;

      if (hyperedge_vertices_in_part_[partOffset] == 0) {
        parts_spanned_[numPartsSpanned++] = vPart;
      }

      ++hyperedge_vertices_in_part_[partOffset];
    }

    // ###
    // update the vertex neighbour structs
    // ###

    for (j = hEdgeOffsets[i]; j < endIndex; ++j) {
      v = pinList[j];
      vertOffset = neighbors_of_vertex_offsets_[v];

      for (ij = 0; ij < numPartsSpanned; ++ij) {
        vPart = parts_spanned_[ij];

        if (neighbors_of_vertex_[vertOffset + vPart] == 0)
          ++number_of_neighboring_parts_[v];

        neighbors_of_vertex_[vertOffset + vPart] = 1;
      }
    }

    k_1cut += ((numPartsSpanned - 1) * hEdgeWeight[i]);
  }

  return k_1cut;
}

void greedy_k_way_refiner::update_adjacent_vertex_stats(int v, int sP,
                                                        int bestDP) {
  int i;
  int j;
  int ij;
  int numVerticesSeen;
  int vertOffset;
  int vert;
  int hEdge;
  int hEdgeOffset;
  int neighOfVOffset;
  int othVOffset;
  int othHedge;
  int reduceNumNeighPartsOn;

  vertOffset = vOffsets[v + 1];
  numVerticesSeen = 0;

  parts_spanned_.resize(number_of_parts_);

  for (j = vOffsets[v]; j < vertOffset; ++j) {

    hEdge = vToHedges[j];
    hEdgeOffset = hEdgeOffsets[hEdge + 1];

    for (ij = hEdgeOffsets[hEdge]; ij < hEdgeOffset; ++ij) {
      vert = pinList[ij];

      if (vertex_seen_[vert] == -1 && vert != v) {
        neighOfVOffset = neighbors_of_vertex_offsets_[vert];

        if (neighbors_of_vertex_[neighOfVOffset + bestDP] == 0)
          ++number_of_neighboring_parts_[vert];

        neighbors_of_vertex_[neighOfVOffset + bestDP] = 1;

        if (partition_vector_[vert] != sP) {
          if (neighbors_of_vertex_[neighOfVOffset + sP] > 0)
            reduceNumNeighPartsOn = 1;
          else
            reduceNumNeighPartsOn = 0;

          neighbors_of_vertex_[neighOfVOffset + sP] = 0;

          othVOffset = vOffsets[vert + 1];

          for (i = vOffsets[vert]; i < othVOffset; ++i) {
            othHedge = vToHedges[i];

            if (hyperedge_vertices_in_part_[hyperedge_vertices_in_part_offsets_[othHedge] + sP] > 0) {
              neighbors_of_vertex_[neighOfVOffset + sP] = 1;
              break;
            }
          }

          if (reduceNumNeighPartsOn && neighbors_of_vertex_[neighOfVOffset + sP] == 0)
            --number_of_neighboring_parts_[vert];
        }

        vertex_seen_[vert] = 1;
        seen_vertices_[numVerticesSeen] = vert;
        ++numVerticesSeen;
      }
    }
  }

  // ###
  // restore the 'seen' vertices structure
  // ###

  for (j = 0; j < numVerticesSeen; ++j)
    vertex_seen_[seen_vertices_[j]] = -1;
}

void greedy_k_way_refiner::refine(serial::hypergraph &h) {
  int totalGain = 0;
  int gain;
  int i;

  load_for_refinement(h);
  build_data_structures();

  if (limit_ < 1.0) {
    set_limit();
  }

  for (i = 0; i < numPartitions; ++i) {
    set_partition_vector(i);
    initialize_data_structures();
    totalGain = 0;

    do {
      gain = greedy_pass();
      totalGain += gain;
    } while (gain > 0);

    partitionCutsizes[i] -= totalGain;
  }

  destroy_data_structures();
}

void greedy_k_way_refiner::rebalance(serial::hypergraph &h) {
  int totalGain = 0;
  int gain;
  int i;

  load_for_refinement(h);
  build_data_structures();

  for (i = 0; i < numPartitions; ++i) {
    set_partition_vector(i);
    initialize_data_structures();

    totalGain = rebalancing_pass();

    do {
      gain = greedy_pass();
      totalGain += gain;
    } while (gain > 0);

    partitionCutsizes[i] -= totalGain;
  }

  destroy_data_structures();

#ifdef DEBUG_REFINER
  h.checkPartitions(numParts, maxPartWt);
#endif
}

int greedy_k_way_refiner::greedy_pass() {
  int i;
  int j;
  int ij;
  int v;
  int sP;
  int gain;
  int vGain;
  int posGain;
  int hEdge;
  int bestMove;
  int randomNum;
  int vertexWt;
  int vertOffset;
  int vNeighOffset;
  int hEdgeOffset;
  int neighOfVOffset;
  int numNonPos = 0;

  double currImbalance = 0;
  double bestImbalance;
  double posImbalance;

  for (i = 0; i < numVertices; ++i) {
    vertices_[i] = i;
  }

  for (i = 0; i < number_of_parts_; ++i) {
    currImbalance += fabs(part_weights_[i] - average_part_weight_);
  }

  i = numVertices;
  gain = 0;

  do {
    randomNum = RANDOM(0, i);
    v = vertices_[randomNum];
    vertexWt = vWeight[v];
    sP = partition_vector_[v];
    vGain = 0;
    bestImbalance = currImbalance;
    bestMove = -1;

    if (number_of_neighboring_parts_[v] > 1) {
      vNeighOffset = neighbors_of_vertex_offsets_[v];

      for (j = 0; j < number_of_parts_; ++j) {
        if (j != sP && neighbors_of_vertex_[vNeighOffset + j] > 0) {
          if (part_weights_[j] + vertexWt <= maximum_part_weight_) {
            posGain = 0;
            vertOffset = vOffsets[v + 1];

            for (ij = vOffsets[v]; ij < vertOffset; ++ij) {
              hEdge = vToHedges[ij];
              hEdgeOffset = hyperedge_vertices_in_part_offsets_[hEdge];

              if (hyperedge_vertices_in_part_[hEdgeOffset + sP] == 1)
                posGain += hEdgeWeight[hEdge];

              if (hyperedge_vertices_in_part_[hEdgeOffset + j] == 0)
                posGain -= hEdgeWeight[hEdge];
            }

            posImbalance =
                currImbalance + fabs(part_weights_[sP] - (vertexWt +
                                                        average_part_weight_));
            posImbalance += fabs((part_weights_[j] + vertexWt) -
                                 average_part_weight_);
            posImbalance -= fabs(part_weights_[sP] - average_part_weight_);
            posImbalance -= fabs(part_weights_[j] - average_part_weight_);

            if ((posGain > vGain) ||
                (posGain == vGain && posImbalance < bestImbalance)) {
              vGain = posGain;
              bestMove = j;
              bestImbalance = posImbalance;
            }
          }
        }
      }

      if (bestMove > -1) {
        vertOffset = vOffsets[v + 1];

        // ###
        // update the moved vertices' stats
        // ###

        neighOfVOffset = neighbors_of_vertex_offsets_[v];

        if (neighbors_of_vertex_[neighOfVOffset + bestMove] == 0)
          ++number_of_neighboring_parts_[v];

        neighbors_of_vertex_[neighOfVOffset + bestMove] = 1;
        neighbors_of_vertex_[neighOfVOffset + sP] = 0;

        for (j = vOffsets[v]; j < vertOffset; ++j) {
          // ###
          // update the hyperedge stats: (vInPart etc.)
          // ###

          hEdge = vToHedges[j];
          hEdgeOffset = hyperedge_vertices_in_part_offsets_[hEdge];

          --hyperedge_vertices_in_part_[hEdgeOffset + sP];
          ++hyperedge_vertices_in_part_[hEdgeOffset + bestMove];

          if (hyperedge_vertices_in_part_[hEdgeOffset + sP] > 0)
            neighbors_of_vertex_[neighOfVOffset + sP] = 1;
        }

        if (neighbors_of_vertex_[neighOfVOffset + sP] == 0)
          --number_of_neighboring_parts_[v];

        // ###
        // update the adj vertices stats:
        // (num neighbours in part etc.)
        // ###

        update_adjacent_vertex_stats(v, sP, bestMove);

        // ###
        // update other structs
        // ###

        partition_vector_[v] = bestMove;
        part_weights_[sP] -= vertexWt;
        part_weights_[bestMove] += vertexWt;
        currImbalance = bestImbalance;

        // ###
        // finally, update the gain...
        // ###

        gain += vGain;

        if (vGain <= 0)
          ++numNonPos;
        else
          numNonPos = 0;

        if (limit_ < 1.0) {
          if (numNonPos > number_of_non_positive_moves_)
            break;
        }
      }
    }

    std::swap(vertices_[randomNum], vertices_[i - 1]);
    --i;
  } while (i > 0);

  return gain;
}

int greedy_k_way_refiner::rebalancing_pass() {
  int i;
  int j;

  int part;
  int numOverWeight = 0;

  dynamic_array<int> overWeight(number_of_parts_);

  for (i = 0; i < number_of_parts_; ++i) {
    if (part_weights_[i] <= maximum_part_weight_)
      overWeight[i] = 0;

    else {
      overWeight[i] = 1;
      ++numOverWeight;
    }
  }

  if (numOverWeight == 0)
    return 0;

  int vertex;
  int hEdge;
  int gain = 0;
  int vGain;
  int posGain;
  int bestMove;
  int vertexWt;
  int vertOffset;
  int hEdgeOffset;
  int neighOfVOffset;

  VNodePtr nodePtr;

  VNodePtrArray verticesInParts(number_of_parts_);
  VNodePtrArray vertexNodes(numVertices);

  for (i = 0; i < numVertices; ++i)
    vertices_[i] = i;

  for (i = 0; i < number_of_parts_; ++i)
    verticesInParts[i] = nullptr;

  Funct::randomPermutation(vertices_.data(), numVertices);

  for (i = 0; i < numVertices; ++i) {
    vertex = vertices_[i];
    part = partition_vector_[vertex];

    if (overWeight[part] == 0) {
      vertexNodes[vertex] = nullptr;
    } else {
      vertexNodes[vertex] = new VNode;
      vertexNodes[vertex]->vertexID = vertex;
      vertexNodes[vertex]->next = verticesInParts[part];
      verticesInParts[part] = vertexNodes[vertex];
    }
  }

  part = find_heaviest_overweight();

  while (part > -1) {
#ifdef DEBUG_REFINER
    assert(part >= 0 && part < numParts);
#endif

    nodePtr = verticesInParts[part];
    verticesInParts[part] = nodePtr->next;
    vertex = nodePtr->vertexID;
    vertexWt = vWeight[vertex];
    vGain = -LARGE_CONSTANT;
    bestMove = -1;

#ifdef DEBUG_REFINER
    assert(vertex >= 0 && vertex < numVertices);
    assert(partitionVector[vertex] == part);
#endif

    for (i = 0; i < number_of_parts_; ++i) {
      if (i != part && overWeight[i] == 0) {
        if (part_weights_[i] + vertexWt <= maximum_part_weight_) {
          posGain = 0;
          vertOffset = vOffsets[vertex + 1];

          for (j = vOffsets[vertex]; j < vertOffset; ++j) {
            hEdge = vToHedges[j];
            hEdgeOffset = hyperedge_vertices_in_part_offsets_[hEdge];

            if (hyperedge_vertices_in_part_[hEdgeOffset + part] == 1)
              posGain += hEdgeWeight[hEdge];

            if (hyperedge_vertices_in_part_[hEdgeOffset + i] == 0)
              posGain -= hEdgeWeight[hEdge];
          }

          if (posGain > vGain) {
            vGain = posGain;
            bestMove = i;
          }
        }
      }
    }
#ifdef DEBUG_REFINER
    assert(bestMove != -1);
#endif

    vertOffset = vOffsets[vertex + 1];

    // ###
    // update the moved vertices' stats
    // ###

    neighOfVOffset = neighbors_of_vertex_offsets_[vertex];

    if (neighbors_of_vertex_[neighOfVOffset + bestMove] == 0)
      ++number_of_neighboring_parts_[vertex];

    neighbors_of_vertex_[neighOfVOffset + bestMove] = 1;
    neighbors_of_vertex_[neighOfVOffset + part] = 0;

    for (i = vOffsets[vertex]; i < vertOffset; ++i) {
      // ###
      // update the hyperedge stats: (vInPart etc.)
      // ###

      hEdge = vToHedges[i];
      hEdgeOffset = hyperedge_vertices_in_part_offsets_[hEdge];

      --hyperedge_vertices_in_part_[hEdgeOffset + part];
      ++hyperedge_vertices_in_part_[hEdgeOffset + bestMove];

      if (hyperedge_vertices_in_part_[hEdgeOffset + part] > 0)
        neighbors_of_vertex_[neighOfVOffset + part] = 1;
    }

    if (neighbors_of_vertex_[neighOfVOffset + part] == 0)
      --number_of_neighboring_parts_[vertex];

    // ###
    // update the adj vertices stats:
    // (num neighbours in part etc.)
    // ###

    update_adjacent_vertex_stats(vertex, part, bestMove);

    // ###
    // update other structs
    // ###

    partition_vector_[vertex] = bestMove;
    part_weights_[part] -= vertexWt;
    part_weights_[bestMove] += vertexWt;

    if (part_weights_[part] <= maximum_part_weight_)
      overWeight[part] = 0;

    // ###
    // finally, update the gain...
    // ###

    gain += vGain;
    part = find_heaviest_overweight();
  }

  return gain;
}

}  // namespace serial
}  // namespace parkway
