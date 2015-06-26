// ### RestrFCCoarsener.cpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 4/1/2005: Last Modified
//
// ###
#include "coarseners/serial/restrictive_first_choice_coarsener.hpp"
#include "utility/logging.hpp"

namespace parkway {
namespace serial {

restrictive_first_choice_coarsener::restrictive_first_choice_coarsener(
    int _min, int _maxwt, double r, int fanOut, int dbWt)
    : restrictive_coarsener(_min, _maxwt, r) {
  util_fan_out_ = fanOut;
  divide_by_cluster_weight_ = dbWt;
}

restrictive_first_choice_coarsener::~restrictive_first_choice_coarsener() {
}

void restrictive_first_choice_coarsener::display_options() const {
  info("[Serial Coarsener: Restrictive First Choice]\n"
       "-- Reduction ratio: %.2f\n"
       "-- Minimum number of nodes: %i\n"
       "-- Util fan out: %i\n\n",
       reduction_ratio_, minimum_nodes_, util_fan_out_);
}

serial::hypergraph *restrictive_first_choice_coarsener::coarsen(const serial::hypergraph &h) {
  loadHypergraphForRestrCoarsening(h);

  // ###
  // back off if the hypergraph too small
  // ###

  if (numVertices < minimum_nodes_)
    return nullptr;

  int numNotMatched = numVertices;
  int bestMatch;
  int bestMatchWt = -1;
  int numNeighbours;
  int limitOnIndex;

  int candVertex;
  int candVertexEntry;

  int i;
  int j;
  int ij;

  int v;
  int vPart;
  int coarseIndex;

  int endVerOffset;
  int endHedgeOffset;
  int hEdge;

  double reducedBy;
  double maxMatchMetric;
  double metricVal;

  dynamic_array<int> neighVerts;
  dynamic_array<int> neighPairWts;
  dynamic_array<int> vertices(numVertices);
  dynamic_array<int> vertexAdjEntry(numVertices);
  dynamic_array<int> coarseWts;
  dynamic_array<int> coarsePVector;

  dynamic_array<double> connectVals;

  bit_field toLoad(numHedges);
  toLoad.set();

  for (i = 0; i < numVertices; ++i) {
    vertices[i] = i;
    vertexAdjEntry[i] = -1;
  }

  if (currPercentile < 100)
    compute_hyperedges_to_load(toLoad);

  Funct::randomPermutation(vertices.data(), numVertices);

  metricVal = static_cast<double>(numVertices) / reduction_ratio_;
  limitOnIndex = numVertices - static_cast<int>(floor(metricVal - 1.0));

  i = 0;
  coarseIndex = 0;

  for (; i < numVertices; ++i) {
    if (matchVector[vertices[i]] == -1) {
      v = vertices[i];
      vPart = partitionVectors[v];
      bestMatch = -1;
      maxMatchMetric = 0;
      numNeighbours = 0;

      endVerOffset = vOffsets[v + 1];

      for (j = vOffsets[v]; j < endVerOffset; ++j) {
        hEdge = vToHedges[j];

        if (toLoad(hEdge) == 1) {
          endHedgeOffset = hEdgeOffsets[hEdge + 1];

          for (ij = hEdgeOffsets[hEdge]; ij < endHedgeOffset; ++ij) {
            candVertex = pinList[ij];

            if (candVertex != v && partitionVectors[candVertex] == vPart) {
              candVertexEntry = vertexAdjEntry[candVertex];

              if (candVertexEntry == -1) {
                neighVerts[numNeighbours] = candVertex;

                if (matchVector[candVertex] == -1)
                  neighPairWts[numNeighbours] = vWeight[v] + vWeight[candVertex];
                else
                  neighPairWts[numNeighbours] = vWeight[v] + coarseWts[matchVector[candVertex]];

                if (util_fan_out_)
                  connectVals[numNeighbours] =
                      static_cast<double>(hEdgeWeight[hEdge]) /
                      (endHedgeOffset - (hEdgeOffsets[hEdge] + 1));
                else
                  connectVals[numNeighbours] =
                      static_cast<double>(hEdgeWeight[hEdge]);

                vertexAdjEntry[candVertex] = numNeighbours;
                ++numNeighbours;
              } else {
                if (util_fan_out_)
                  connectVals[candVertexEntry] +=
                      static_cast<double>(hEdgeWeight[hEdge]) /
                      (endHedgeOffset - (hEdgeOffsets[hEdge] + 1));
                else
                  connectVals[candVertexEntry] +=
                      static_cast<double>(hEdgeWeight[hEdge]);
              }
            }
          }
        }
      }

      if (numNeighbours == 0) {
        // match vertex v as singleton as not connected

        matchVector[v] = coarseIndex;
        coarsePVector[coarseIndex] = vPart;
        coarseWts[coarseIndex++] = vWeight[v];
        --numNotMatched;
      } else {
        // else need to pick the best match

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

          matchVector[v] = coarseIndex;
          coarsePVector[coarseIndex] = vPart;
          coarseWts[coarseIndex++] = vWeight[v];
          --numNotMatched;
        } else {
          if (matchVector[bestMatch] == -1) {
            // match with another unmatched

            matchVector[v] = coarseIndex;
            matchVector[bestMatch] = coarseIndex;
            coarsePVector[coarseIndex] = vPart;
            coarseWts[coarseIndex++] = bestMatchWt;
            numNotMatched -= 2;
          } else {
            // match with existing cluster of coarse vertices

            matchVector[v] = matchVector[bestMatch];
            coarseWts[matchVector[v]] += vWeight[v];
            --numNotMatched;
          }
        }
      }

      // check if hypergraph sufficiently shrunk

      if (i > limitOnIndex) {
        reducedBy = static_cast<double>(numVertices) /
            (numNotMatched + coarseIndex);

        if (reducedBy > reduction_ratio_)
          break;
      }
    }
  }

  if (i == numVertices) {
    // cannot sufficiently reduce hypergraph
    // back off
    return nullptr;
  }

  // match remaining unmatched vertices as singletons

  for (; i < numVertices; ++i) {
    v = vertices[i];

    if (matchVector[v] == -1) {
      matchVector[v] = coarseIndex;
      coarsePVector[coarseIndex] = partitionVectors[v];
      coarseWts[coarseIndex++] = vWeight[v];
    }
  }

  coarseWts.resize(coarseIndex);
  coarsePVector.resize(coarseIndex);

  return (build_coarse_hypergraph(coarseWts, coarsePVector, coarseIndex,
                                  h.total_weight()));
}

}  // namespace serial
}  // namespace parkway
