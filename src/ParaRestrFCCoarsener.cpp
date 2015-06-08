
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
                                           int divByLen, ostream &out)
    : ParaRestrCoarsener(rank, nProcs, nParts, out) {
  vertexVisitOrder = verVisOrder;
  divByCluWt = divByWt;
  divByHedgeLen = divByLen;
  limitOnIndexDuringCoarsening = 0;
}

ParaRestrFCCoarsener::~ParaRestrFCCoarsener() {}

void ParaRestrFCCoarsener::dispCoarseningOptions() const {
  switch (dispOption) {
  case SILENT:

    break;

  default:

    out_stream << "|--- PARA_RESTR_C:" << endl
               << "|- PFC:"
               << " r = " << reductionRatio << " min = " << minNodes
               << " vvo = ";
    printVisitOrder(vertexVisitOrder);
    out_stream << " divWt = " << divByCluWt << " divLen = " << divByHedgeLen
               << endl
               << "|" << endl;
    break;
  }
}

void ParaRestrFCCoarsener::buildAuxiliaryStructs(int numPins, double aveVertDeg,
                                                 double aveHedgeSize) {}

void ParaRestrFCCoarsener::releaseMemory() {
  hEdgeWeight.setLength(0);
  hEdgeOffset.setLength(0);
  locPinList.setLength(0);

  vToHedgesOffset.setLength(0);
  vToHedgesList.setLength(0);
  allocHedges.setLength(0);

  freeMemory();
}

ParaHypergraph *ParaRestrFCCoarsener::coarsen(ParaHypergraph &h,
                                              MPI_Comm comm) {
  loadHyperGraph(h, comm);

  if (totalVertices < minNodes || h.dontCoarsen()) {
    return NULL;
  }

  int i;
  int j;

  int index = 0;
  int numNotMatched = numLocalVertices;
  int aveVertexWt = static_cast<int>(
      ceil(static_cast<double>(totalHypergraphWt) / totalVertices));
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

  DynamicArray<int> neighVerts;
  DynamicArray<int> neighPairWts;
  DynamicArray<double> connectVals;
  DynamicArray<int> vertexAdjEntry(numLocalVertices);
  DynamicArray<int> vertices(numLocalVertices);

  permuteVerticesArray(vertices.getArray(), numLocalVertices);

  for (i = 0; i < numLocalVertices; ++i)
    vertexAdjEntry[i] = -1;

  if (dispOption > 1) {
    for (i = 0; i < numLocalVertices; ++i) {
      if (vWeight[i] > maxLocWt)
        maxLocWt = vWeight[i];
    }

    MPI_Reduce(&maxLocWt, &maxWt, 1, MPI_INT, MPI_MAX, 0, comm);

    if (myRank == 0) {
      out_stream << " " << maxVertexWt << " " << maxWt << " " << aveVertexWt
                 << " ";
    }
  }

  metricVal = static_cast<double>(numLocalVertices) / reductionRatio;

#ifdef DEBUG_COARSENER
  assert(clusterWeights == NULL);
#endif

  clusterWeights = new DynamicArray<int>(1024);
  pVector = new DynamicArray<int>(1024);

  limitOnIndexDuringCoarsening =
      numLocalVertices - static_cast<int>(floor(metricVal - 1.0));
  clusterIndex = 0;

  for (; index < numLocalVertices; ++index) {
    if (matchVector[vertices[index]] == -1) {
      vertex = vertices[index];
#ifdef DEBUG_COARSENER
      assert(vertex >= 0 && vertex < numLocalVertices);
#endif
      vertexPart = partitionVector[vertex];
      endOffset1 = vToHedgesOffset[vertex + 1];
      numNeighbours = 0;

      for (i = vToHedgesOffset[vertex]; i < endOffset1; ++i) {
        hEdge = vToHedgesList[i];
        hEdgeWt = hEdgeWeight[hEdge];
        endOffset2 = hEdgeOffset[hEdge + 1];
        hEdgeLen = endOffset2 - hEdgeOffset[hEdge];

        for (j = hEdgeOffset[hEdge]; j < endOffset2; ++j) {
#ifdef DEBUG_COARSENER
          assert(locPinList[j] >= 0 && locPinList[j] < numLocalVertices);
#endif
          neighVertex = locPinList[j];

          if (neighVertex != vertex &&
              partitionVector[neighVertex] == vertexPart) {
            // ###
            // now try not to check the weight before checking the vertex
            // ###

            neighVertexEntry = vertexAdjEntry[neighVertex];

            if (neighVertexEntry == -1) {
              if (matchVector[neighVertex] == -1)
                neighPairWts.assign(numNeighbours,
                                    vWeight[vertex] + vWeight[neighVertex]);
              else
                neighPairWts.assign(
                    numNeighbours,
                    vWeight[vertex] +
                        (*clusterWeights)[matchVector[neighVertex]]);

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

        matchVector[vertex] = clusterIndex;
        pVector->assign(clusterIndex, partitionVector[vertex]);
        clusterWeights->assign(clusterIndex++, vWeight[vertex]);
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

          matchVector[vertex] = clusterIndex;
          pVector->assign(clusterIndex, partitionVector[vertex]);
          clusterWeights->assign(clusterIndex++, vWeight[vertex]);
          --numNotMatched;
        } else {
#ifdef DEBUG_COARSENER
          assert(bestMatch >= 0 && bestMatch < numLocalVertices);
#endif
          if (matchVector[bestMatch] == -1) {
            // match with another unmatched

            matchVector[vertex] = clusterIndex;
            matchVector[bestMatch] = clusterIndex;
            pVector->assign(clusterIndex, partitionVector[vertex]);
            clusterWeights->assign(clusterIndex++, bestMatchWt);
            numNotMatched -= 2;
          } else {
            // match with existing cluster of coarse vertices

            matchVector[vertex] = matchVector[bestMatch];
            (*clusterWeights)[matchVector[vertex]] += vWeight[vertex];
            --numNotMatched;
          }
        }
      }

      // check if hypergraph sufficiently shrunk

      if (index > limitOnIndexDuringCoarsening) {
        reducedBy = static_cast<double>(numLocalVertices) /
                    (numNotMatched + clusterIndex);

        if (reducedBy > reductionRatio)
          break;
      }
    }
  }

  for (; index < numLocalVertices; ++index) {
    vertex = vertices[index];

    if (matchVector[vertex] == -1) {
      matchVector[vertex] = clusterIndex;
      pVector->assign(clusterIndex, partitionVector[vertex]);
      clusterWeights->assign(clusterIndex++, vWeight[vertex]);
    }
  }

#ifdef DEBUG_COARSENER
  for (i = 0; i < numLocalVertices; ++i) {
    assert(matchVector[i] != -1);
    assert(matchVector[i] >= 0 && matchVector[i] < clusterIndex);
  }
#endif

  setClusterIndices(comm);

  if (static_cast<double>(totalVertices) / totalClusters <
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
    Funct::qsortByAnotherArray(0, nLocVerts - 1, verts, vWeight, INC);
    break;

  case DECREASING_WEIGHT_ORDER:
    for (i = 0; i < nLocVerts; ++i) {
      verts[i] = i;
    }
    Funct::qsortByAnotherArray(0, nLocVerts - 1, verts, vWeight, DEC);
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

  if (myRank == numProcs - 1)
    totalClusters = myMinCluIndex;

  myMinCluIndex -= clusterIndex;

  MPI_Bcast(&totalClusters, 1, MPI_INT, numProcs - 1, comm);

  for (i = 0; i < numLocalVertices; ++i) {
#ifdef DEBUG_COARSENER
    assert(matchVector[i] != -1);
#endif
    matchVector[i] += myMinCluIndex;
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
