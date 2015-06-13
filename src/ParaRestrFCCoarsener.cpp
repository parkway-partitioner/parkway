#ifndef _PARA_RESTR_FCC_CPP
#define _PARA_RESTR_FCC_CPP

// ### ParaRestrFCCoarsener.cpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 31/12/2004: Last Modified
//
// ###

#include "ParaRestrFCCoarsener.hpp"

ParaRestrFCCoarsener::ParaRestrFCCoarsener(int rank, int nProcs, int nParts,
                                           int verVisOrder, int divByWt,
                                           int divByLen, std::ostream &out)
    : ParaRestrCoarsener(rank, nProcs, nParts, out) {
  vertexVisitOrder = verVisOrder;
  divByCluWt = divByWt;
  divByHedgeLen = divByLen;
  limitOnIndexDuringCoarsening = 0;
}

ParaRestrFCCoarsener::~ParaRestrFCCoarsener() {}

void ParaRestrFCCoarsener::dispCoarseningOptions() const {
  switch (display_options_) {
  case SILENT:

    break;

  default:

    out_stream << "|--- PARA_RESTR_C:" << std::endl
               << "|- PFC:"
               << " r = " << reductionRatio << " min = " << minNodes
               << " vvo = ";
    printVisitOrder(vertexVisitOrder);
    out_stream << " divWt = " << divByCluWt << " divLen = " << divByHedgeLen
               << std::endl
               << "|" << std::endl;
    break;
  }
}

void ParaRestrFCCoarsener::buildAuxiliaryStructs(int numPins, double aveVertDeg,
                                                 double aveHedgeSize) {}

void ParaRestrFCCoarsener::release_memory() {
  hyperedge_weights_.reserve(0);
  hyperedge_offsets_.reserve(0);
  local_pin_list_.reserve(0);

  vertex_to_hyperedges_offset_.reserve(0);
  vertex_to_hyperedges_.reserve(0);
  allocated_hyperedges_.reserve(0);

  free_memory();
}

hypergraph *ParaRestrFCCoarsener::coarsen(hypergraph &h, MPI_Comm comm) {
  load(h, comm);

  if (number_of_vertices_ < minNodes || h.dont_coarsen()) {
    return nullptr;
  }

  int i;
  int j;

  int index = 0;
  int numNotMatched = number_of_local_vertices_;
  int aveVertexWt = static_cast<int>(
      ceil(static_cast<double>(totalHypergraphWt) / number_of_vertices_));
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

  dynamic_array<int> neighVerts;
  dynamic_array<int> neighPairWts;
  dynamic_array<double> connectVals;
  dynamic_array<int> vertexAdjEntry(number_of_local_vertices_);
  dynamic_array<int> vertices(number_of_local_vertices_);

  permuteVerticesArray(vertices.data(), number_of_local_vertices_);

  for (i = 0; i < number_of_local_vertices_; ++i)
    vertexAdjEntry[i] = -1;

  if (display_options_ > 1) {
    for (i = 0; i < number_of_local_vertices_; ++i) {
      if (vertex_weights_[i] > maxLocWt)
        maxLocWt = vertex_weights_[i];
    }

    MPI_Reduce(&maxLocWt, &maxWt, 1, MPI_INT, MPI_MAX, 0, comm);

    if (rank_ == 0) {
      out_stream << " " << maxVertexWt << " " << maxWt << " " << aveVertexWt
                 << " ";
    }
  }

  metricVal = static_cast<double>(number_of_local_vertices_) / reductionRatio;

#ifdef DEBUG_COARSENER
  assert(clusterWeights == nullptr);
#endif

  clusterWeights = new dynamic_array<int>(1024);
  pVector = new dynamic_array<int>(1024);

  limitOnIndexDuringCoarsening =
      number_of_local_vertices_ - static_cast<int>(floor(metricVal - 1.0));
  clusterIndex = 0;

  for (; index < number_of_local_vertices_; ++index) {
    if (match_vector_[vertices[index]] == -1) {
      vertex = vertices[index];
#ifdef DEBUG_COARSENER
      assert(vertex >= 0 && vertex < numLocalVertices);
#endif
      vertexPart = partitionVector[vertex];
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
              partitionVector[neighVertex] == vertexPart) {
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
                        (*clusterWeights)[match_vector_[neighVertex]]);

              neighVerts.assign(numNeighbours, neighVertex);

              if (divByHedgeLen)
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
              if (divByHedgeLen)
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

        match_vector_[vertex] = clusterIndex;
        pVector->assign(clusterIndex, partitionVector[vertex]);
        clusterWeights->assign(clusterIndex++, vertex_weights_[vertex]);
        --numNotMatched;
      } else {
        // ###
        // pick best match
        // ###

        bestMatch = -1;
        maxMatchMetric = 0.0;

        for (j = 0; j < numNeighbours; ++j) {
          if (neighPairWts[j] < maxVertexWt) {
            metricVal = connectVals[j];

            if (divByCluWt)
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

          match_vector_[vertex] = clusterIndex;
          pVector->assign(clusterIndex, partitionVector[vertex]);
          clusterWeights->assign(clusterIndex++, vertex_weights_[vertex]);
          --numNotMatched;
        } else {
#ifdef DEBUG_COARSENER
          assert(bestMatch >= 0 && bestMatch < numLocalVertices);
#endif
          if (match_vector_[bestMatch] == -1) {
            // match with another unmatched

            match_vector_[vertex] = clusterIndex;
            match_vector_[bestMatch] = clusterIndex;
            pVector->assign(clusterIndex, partitionVector[vertex]);
            clusterWeights->assign(clusterIndex++, bestMatchWt);
            numNotMatched -= 2;
          } else {
            // match with existing cluster of coarse vertices

            match_vector_[vertex] = match_vector_[bestMatch];
            (*clusterWeights)[match_vector_[vertex]] += vertex_weights_[vertex];
            --numNotMatched;
          }
        }
      }

      // check if hypergraph sufficiently shrunk

      if (index > limitOnIndexDuringCoarsening) {
        reducedBy = static_cast<double>(number_of_local_vertices_) /
                    (numNotMatched + clusterIndex);

        if (reducedBy > reductionRatio)
          break;
      }
    }
  }

  for (; index < number_of_local_vertices_; ++index) {
    vertex = vertices[index];

    if (match_vector_[vertex] == -1) {
      match_vector_[vertex] = clusterIndex;
      pVector->assign(clusterIndex, partitionVector[vertex]);
      clusterWeights->assign(clusterIndex++, vertex_weights_[vertex]);
    }
  }

#ifdef DEBUG_COARSENER
  for (i = 0; i < numLocalVertices; ++i) {
    assert(matchVector[i] != -1);
    assert(matchVector[i] >= 0 && matchVector[i] < clusterIndex);
  }
#endif

  setClusterIndices(comm);

  if (static_cast<double>(number_of_vertices_) / totalClusters <
      MIN_ALLOWED_REDUCTION_RATIO)
    stopCoarsening = 1;
  else
    stopCoarsening = 0;

  // ###
  // now construct the coarse hypergraph using the matching vector
  // ###

  return (contractHyperedges(h, comm));
}

void ParaRestrFCCoarsener::permuteVerticesArray(int *verts, int nLocVerts) {
  int i;

  switch (vertexVisitOrder) {
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
    Funct::randomPermutation(verts, nLocVerts);
    break;

  case INCREASING_WEIGHT_ORDER:
    for (i = 0; i < nLocVerts; ++i) {
      verts[i] = i;
    }
    Funct::qsortByAnotherArray(0, nLocVerts - 1, verts, vertex_weights_, INC);
    break;

  case DECREASING_WEIGHT_ORDER:
    for (i = 0; i < nLocVerts; ++i) {
      verts[i] = i;
    }
    Funct::qsortByAnotherArray(0, nLocVerts - 1, verts, vertex_weights_, DEC);
    break;

  default:
    for (i = 0; i < nLocVerts; ++i) {
      verts[i] = i;
    }
    Funct::randomPermutation(verts, nLocVerts);
    break;
  }
}

void ParaRestrFCCoarsener::setClusterIndices(MPI_Comm comm) {
  int i;

  MPI_Scan(&clusterIndex, &myMinCluIndex, 1, MPI_INT, MPI_SUM, comm);

  if (rank_ == processors_ - 1)
    totalClusters = myMinCluIndex;

  myMinCluIndex -= clusterIndex;

  MPI_Bcast(&totalClusters, 1, MPI_INT, processors_ - 1, comm);

  for (i = 0; i < number_of_local_vertices_; ++i) {
#ifdef DEBUG_COARSENER
    assert(matchVector[i] != -1);
#endif
    match_vector_[i] += myMinCluIndex;
  }
}

void ParaRestrFCCoarsener::printVisitOrder(int variable) const {
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

#endif
