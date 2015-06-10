
#ifndef _PARA_FCCOARSENER_CPP
#define _PARA_FCCOARSENER_CPP

// ### ParaFCCoarsener.cpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 3/2/2005: Last Modified
//
// ###

#include "ParaFCCoarsener.hpp"
#include <iostream>
#include "data_structures/internal/table_utils.hpp"
#include "data_structures/match_request_table.hpp"
#include "data_structures/map_to_pos_int.hpp"

ParaFCCoarsener::ParaFCCoarsener(int rank, int nProcs, int nParts,
                                 int vertVisOrder, int matchReqOrder,
                                 int divByWt, int divByLen, std::ostream &out)
    : ParaCoarsener(rank, nProcs, nParts, out) {
  vertexVisitOrder = vertVisOrder;
  matchRequestVisitOrder = matchReqOrder;
  divByCluWt = divByWt;
  divByHedgeLen = divByLen;
  limitOnIndexDuringCoarsening = 0;

  table = nullptr;
}

ParaFCCoarsener::~ParaFCCoarsener() {
  DynaMem::deletePtr<ds::match_request_table>(table);
}

void ParaFCCoarsener::dispCoarseningOptions() const {
  switch (dispOption) {
  case SILENT:

    break;

  default:

    out_stream << "|--- PARA_C:" << std::endl
               << "|- PFC:"
               << " r = " << reductionRatio << " min = " << minNodes
               << " vvo = ";
    printVisitOrder(vertexVisitOrder);
    out_stream << " mvo = ";
    printVisitOrder(matchRequestVisitOrder);
    out_stream << " divWt = " << divByCluWt << " divLen = " << divByHedgeLen
               << std::endl
               << "|" << std::endl;
    break;
  }
}

void ParaFCCoarsener::buildAuxiliaryStructs(int numTotPins, double aveVertDeg,
                                            double aveHedgeSize) {
  // ###
  // build the ds::match_request_table
  // ###

  int i =
      static_cast<int>(ceil(static_cast<double>(numTotPins) / aveVertDeg));

  table = new ds::match_request_table(ds::internal::table_utils::table_size(i / numProcs));
}

void ParaFCCoarsener::releaseMemory() {
  hEdgeWeight.reserve(0);
  hEdgeOffset.reserve(0);
  locPinList.reserve(0);

  vToHedgesOffset.reserve(0);
  vToHedgesList.reserve(0);
  allocHedges.reserve(0);
  clusterWeights.reserve(0);

  freeMemory();
}

ParaHypergraph *ParaFCCoarsener::coarsen(ParaHypergraph &h, MPI_Comm comm) {
  loadHyperGraph(h, comm);

#ifdef MEM_CHECK
  MPI_Barrier(comm);
  write_log(myRank, "[begin PFCC]: usage: %f", MemoryTracker::usage());
  Funct::printMemUse(myRank, "[begin PFCC]");
#endif

  if (totalVertices < minNodes || h.dontCoarsen()) {
    return nullptr;
  }

  int i;
  int j;

  int index = 0;
  int numNotMatched = numLocalVertices;
  int maxLocWt = 0;
  int maxWt;
  int aveVertexWt = static_cast<int>(
      ceil(static_cast<double>(totalHypergraphWt) / totalVertices));
  int bestMatch;
  int bestMatchWt = -1;
  int numVisited;
  int vPerProc = totalVertices / numProcs;
  int endOffset1;
  int endOffset2;
  int cluWeight;
  int globalVertexIndex;
  int candVwt;
  int nonLocV;
  int candidatV;
  int hEdgeLen;
  int pairWt;
  int neighbourLoc;

  int hEdge;
  int vertex;

  double metric;
  double maxMatchMetric;
  double reducedBy;

  ds::map_to_pos_int matchInfoLoc;

  if (numLocPins < totalVertices / 2)
    matchInfoLoc.create(numLocPins, 1);
  else
    matchInfoLoc.create(totalVertices, 0);

  dynamic_array<int> neighVerts;
  dynamic_array<int> neighCluWts;
  dynamic_array<double> connectVals;
  dynamic_array<int> vertices(numLocalVertices);

  permuteVerticesArray(vertices.data(), numLocalVertices);

  if (dispOption > 1) {
    for (i = 0; i < numLocalVertices; ++i) {
      if (vWeight[i] > maxLocWt)
        maxLocWt = vWeight[i];
    }

    MPI_Reduce(&maxLocWt, &maxWt, 1, MPI_INT, MPI_MAX, 0, comm);

    if (myRank == 0) {
      out_stream << " " << maxVertexWt << " " << maxWt << " " << aveVertexWt
                 << " ";
      out_stream.flush();
    }
  }

  metric = static_cast<double>(numLocalVertices) / reductionRatio;
  limitOnIndexDuringCoarsening =
      numLocalVertices - static_cast<int>(floor(metric - 1.0));
  clusterIndex = 0;
  stopCoarsening = 0;

  for (; index < numLocalVertices; ++index) {
    if (matchVector[vertices[index]] == -1) {
      vertex = vertices[index];
#ifdef DEBUG_COARSENER
      assert(vertex >= 0 && vertex < numLocalVertices);
#endif
      globalVertexIndex = vertex + minVertexIndex;
      endOffset1 = vToHedgesOffset[vertex + 1];
      bestMatch = -1;
      maxMatchMetric = 0.0;
      numVisited = 0;

      for (i = vToHedgesOffset[vertex]; i < endOffset1; ++i) {
        hEdge = vToHedgesList[i];
        endOffset2 = hEdgeOffset[hEdge + 1];
        hEdgeLen = endOffset2 - hEdgeOffset[hEdge];

        for (j = hEdgeOffset[hEdge]; j < endOffset2; ++j) {

          candidatV = locPinList[j];
#ifdef DEBUG_COARSENER
          assert(candidatV >= 0 && candidatV < totalVertices);
#endif
          if (candidatV != globalVertexIndex) {
            /* now try not to check the weight before checking the vertex */

            neighbourLoc = matchInfoLoc.get_careful(candidatV);

            if (neighbourLoc >= 0) {
              if (divByHedgeLen)
                connectVals[neighbourLoc] +=
                    (static_cast<double>(hEdgeWeight[hEdge]) / (hEdgeLen - 1));
              else
                connectVals[neighbourLoc] +=
                    (static_cast<double>(hEdgeWeight[hEdge]));
            } else {
              // ###
              // here compute the cluster weight
              // ###

              if (candidatV >= minVertexIndex && candidatV < maxVertexIndex) {
                // ###
                // candidatV is a local vertex
                // ###

                if (matchVector[candidatV - minVertexIndex] == -1)
                  cluWeight =
                      vWeight[vertex] + vWeight[candidatV - minVertexIndex];
                else if (matchVector[candidatV - minVertexIndex] >=
                         NON_LOCAL_MATCH) {
                  nonLocV =
                      matchVector[candidatV - minVertexIndex] - NON_LOCAL_MATCH;
                  cluWeight = vWeight[vertex] +
                              table->cluster_weight(nonLocV) + aveVertexWt;
                } else
                  cluWeight =
                      vWeight[vertex] +
                      clusterWeights[matchVector[candidatV - minVertexIndex]];
              } else {
                // ###
                // candidatV is not a local vertex or it is a local
                // vertex matched with a non-local one
                // ###

                candVwt = table->cluster_weight(candidatV);

                if (candVwt != -1)
                  cluWeight = vWeight[vertex] + candVwt + aveVertexWt;
                else
                  cluWeight = vWeight[vertex] + aveVertexWt;
              }

#ifdef DEBUG_COARSENER
              assert(numVisited >= 0);
#endif
              if (matchInfoLoc.insert(candidatV, numVisited)) {
                write_log(myRank, "numEntries %d",
                          matchInfoLoc.size());
                write_log(myRank, "using hash %d",
                          matchInfoLoc.use_hash());
                write_log(myRank, "numSlots %d", matchInfoLoc.capacity());
                write_log(myRank, "candidatV %d", candidatV);
                write_log(myRank, "neighbourLoc %d", neighbourLoc);
                assert(0);
              }

              neighVerts.assign(numVisited, candidatV);
              neighCluWts.assign(numVisited, cluWeight);

              if (divByHedgeLen)
                connectVals.assign(numVisited++,
                                   static_cast<double>(hEdgeWeight[hEdge]) /
                                       (hEdgeLen - 1));
              else
                connectVals.assign(numVisited++,
                                   static_cast<double>(hEdgeWeight[hEdge]));
            }
          }
        }
      }

      // ###
      // now choose best match from adjacent vertices
      // visited above
      // ###

      for (i = 0; i < numVisited; ++i) {
        pairWt = neighCluWts[i];
        candidatV = neighVerts[i];

        if (pairWt <= maxVertexWt) {
          metric = connectVals[i];

          if (divByCluWt)
            metric /= pairWt;

          if (metric > maxMatchMetric) {
            maxMatchMetric = metric;
            bestMatch = candidatV;
            bestMatchWt = pairWt;

#ifdef DEBUG_COARSENER
            assert(bestMatch >= 0);
#endif
          }
        }
      }

      matchInfoLoc.clear();
      numVisited = 0;

      if (bestMatch == -1) {
        // ###
        // match as singleton
        // ###

        matchVector[vertex] = clusterIndex;
        clusterWeights.assign(clusterIndex++, vWeight[vertex]);
        --numNotMatched;
      } else {
#ifdef DEBUG_COARSENER
        assert(bestMatch >= 0);
#endif

        if (bestMatch >= minVertexIndex && bestMatch < maxVertexIndex) {
          // ###
          // best match is a local vertex
          // ###

          if (matchVector[bestMatch - minVertexIndex] == -1) {
            matchVector[bestMatch - minVertexIndex] = clusterIndex;
            matchVector[vertex] = clusterIndex;
            clusterWeights.assign(clusterIndex++, bestMatchWt);
            numNotMatched -= 2;
          } else {
            if (matchVector[bestMatch - minVertexIndex] >= NON_LOCAL_MATCH) {
              nonLocV =
                  matchVector[bestMatch - minVertexIndex] - NON_LOCAL_MATCH;
              table->add_local(nonLocV, vertex + minVertexIndex,
                               vWeight[vertex],
                               std::min(nonLocV / vPerProc, numProcs - 1));
#ifdef DEBUG_COARSENER
              assert(std::min(nonLocV / vPerProc, numProcs - 1) != myRank);
#endif
              matchVector[vertex] = NON_LOCAL_MATCH + nonLocV;
              --numNotMatched;
            } else {
              matchVector[vertex] = matchVector[bestMatch - minVertexIndex];
              clusterWeights[matchVector[vertex]] +=
                  vWeight[vertex]; // bestMatchWt;
              --numNotMatched;
            }
          }
        } else {
          // ###
          // best match is not a local vertex
          // ###

          table->add_local(bestMatch, vertex + minVertexIndex, vWeight[vertex],
                           std::min(bestMatch / vPerProc, numProcs - 1));
#ifdef DEBUG_COARSENER
          assert(std::min(bestMatch / vPerProc, numProcs - 1) != myRank);
#endif
          matchVector[vertex] = NON_LOCAL_MATCH + bestMatch;
          --numNotMatched;
        }
      }

      // ###
      // check the amount that the hypergraph has shrunk by
      // ###

      if (index > limitOnIndexDuringCoarsening) {
        reducedBy = static_cast<double>(numLocalVertices) /
                    (numNotMatched + clusterIndex + table->size());
        break;
      }
    }
  }

  matchInfoLoc.destroy();

  // ###
  // now carry over all the unmatched vertices as singletons
  // ###

  for (; index < numLocalVertices; ++index) {
    vertex = vertices[index];

    if (matchVector[vertex] == -1) {
      matchVector[vertex] = clusterIndex;
      clusterWeights.assign(clusterIndex++, vWeight[vertex]);
    }
  }

  // ###
  // now need to prepare and send out the matching requests
  // ###

  for (i = 0; i < 2; ++i) {
    setRequestArrays(i);
    sendFromDataOutArrays(comm); // actually sending requests
    setReplyArrays(i, maxVertexWt);
    sendFromDataOutArrays(comm); // actually sending replies
    processReqReplies();
  }

  setClusterIndices(comm);

  if (static_cast<double>(totalVertices) / totalClusters <
      MIN_ALLOWED_REDUCTION_RATIO) {
    stopCoarsening = 1;
  }

  table->clear();

  // ###
  // now construct the coarse hypergraph using the matching vector
  // ###

  /*
  if(myRank == 0) {
    std::cout << "about to contract hyperedges" << std::endl;
  }
  MPI_Barrier(comm);
  */
  return (contractHyperedges(h, comm));
}

void ParaFCCoarsener::setRequestArrays(int highToLow) {
  int numRequests = table->size();
  int nonLocVertex;
  int cluWt;

  int i;
  int procRank;

  ds::match_request_table::entry *entry_;
  ds::match_request_table::entry **entryArray = table->get_entries();

  for (i = 0; i < numProcs; ++i)
    sendLens[i] = 0;

  for (i = 0; i < numRequests; ++i) {
    entry_ = entryArray[i];

#ifdef DEBUG_COARSENER
    assert(entry_);
#endif

    nonLocVertex = entry_->non_local_vertex();
    cluWt = entry_->cluster_weight();
    procRank = entry_->non_local_process();

#ifdef DEBUG_COARSENER
    assert(procRank < numProcs);
#endif

    if ((cluWt >= 0) && ((highToLow && procRank < myRank) ||
                         (!highToLow && procRank > myRank))) {
      dataOutSets[procRank]->assign(sendLens[procRank]++, nonLocVertex);
      dataOutSets[procRank]->assign(sendLens[procRank]++, cluWt);
    }
  }
}

void ParaFCCoarsener::setReplyArrays(int highToLow, int maxVWt) {
  int j;
  int i;
  int l;

  dynamic_array<int> visitOrder;

  for (i = 0; i < numProcs; ++i)
    sendLens[i] = 0;

  int startOffset = 0;
  int vLocReq;
  int reqCluWt;
  int matchIndex;
  int visitOrderLen;

  for (i = 0; i < numProcs; ++i) {
#ifdef DEBUG_COARSENER
    assert(And(recvLens[i], 0x1) == 0);
#endif

    if (matchRequestVisitOrder == RANDOM_ORDER) {
      visitOrderLen = Shiftr(recvLens[i], 1);
      visitOrder.reserve(visitOrderLen);

      for (j = 0; j < visitOrderLen; ++j)
        visitOrder[j] = j;

      // ###
      // processing match requests in random order
      // ###

      Funct::randomPermutation(visitOrder.data(), visitOrderLen);

      for (l = 0; l < visitOrderLen; ++l) {

        j = Shiftl(visitOrder[l], 1);

        vLocReq = receiveArray[startOffset + j];
        reqCluWt = receiveArray[startOffset + j + 1];

        if (accept(vLocReq, reqCluWt, highToLow, maxVWt)) {
          // ###
          // cross-processor match accepted. inform vertices
          // of their cluster index and cluster weight
          // ###

          matchIndex = matchVector[vLocReq - minVertexIndex];
          dataOutSets[i]->assign(sendLens[i]++, vLocReq);
          dataOutSets[i]->assign(sendLens[i]++, matchIndex);
          dataOutSets[i]->assign(sendLens[i]++, clusterWeights[matchIndex]);
        } else {
          // ###
          // cross-processor match rejected, inform vertices
          // that match rejected
          // ###

          dataOutSets[i]->assign(sendLens[i]++, vLocReq);
          dataOutSets[i]->assign(sendLens[i]++, NO_MATCH);
        }
      }
    } else {
      // ###
      // processing match requests as they arrive
      // ###

      visitOrderLen = recvLens[i];

      for (l = 0; l < visitOrderLen;) {
        vLocReq = receiveArray[startOffset + (l++)];
        reqCluWt = receiveArray[startOffset + (l++)];

        if (accept(vLocReq, reqCluWt, highToLow, maxVWt)) {
          // ###
          // cross-processor match accepted. inform vertices
          // of their cluster index and cluster weight
          // ###

          matchIndex = matchVector[vLocReq - minVertexIndex];
          dataOutSets[i]->assign(sendLens[i]++, vLocReq);
          dataOutSets[i]->assign(sendLens[i]++, matchIndex);
          dataOutSets[i]->assign(sendLens[i]++, clusterWeights[matchIndex]);
        } else {
          // ###
          // cross-processor match rejected, inform vertices
          // that match rejected
          // ###

          dataOutSets[i]->assign(sendLens[i]++, vLocReq);
          dataOutSets[i]->assign(sendLens[i]++, NO_MATCH);
        }
      }
    }
    startOffset += recvLens[i];
  }
}

void ParaFCCoarsener::processReqReplies() {
  int i;
  int j;
  int index;

  int startOffset = 0;
  int vNonLocReq;
  int cluWt;
  int matchIndex;
  int numLocals;
  int *locals;

  ds::match_request_table::entry *entry_;

  for (i = 0; i < numProcs; ++i) {
    j = 0;
    while (j < recvLens[i]) {
      vNonLocReq = receiveArray[startOffset + (j++)];
      matchIndex = receiveArray[startOffset + (j++)];

      if (matchIndex != NO_MATCH) {
        // ###
        // match successful - set the clusterIndex
        // ###

        cluWt = receiveArray[startOffset + (j++)];
        entry_ = table->get_entry(vNonLocReq);

#ifdef DEBUG_COARSENER
        assert(entry_);
#endif

        entry_->set_cluster_index(matchIndex);
        entry_->set_cluster_weight(cluWt);
      } else {
        // ###
        // match not successful - match requesting
        // vertices into a cluster
        // ###

        entry_ = table->get_entry(vNonLocReq);
        locals = entry_->local_vertices_array();
        numLocals = entry_->number_local();
        entry_->set_cluster_index(MATCHED_LOCALLY);

        for (index = 0; index < numLocals; ++index)
          matchVector[locals[index] - minVertexIndex] = clusterIndex;

        clusterWeights.assign(clusterIndex++, entry_->cluster_weight());
      }
    }
    startOffset += recvLens[i];
  }
}

void ParaFCCoarsener::permuteVerticesArray(int *verts, int nLocVerts) {
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

void ParaFCCoarsener::setClusterIndices(MPI_Comm comm) {
  dynamic_array<int> numClusters(numProcs);
  dynamic_array<int> startIndex(numProcs);

  MPI_Allgather(&clusterIndex, 1, MPI_INT, numClusters.data(), 1, MPI_INT,
                comm);

  int index = 0;
  int i;

  ds::match_request_table::entry *entry_;
  ds::match_request_table::entry **entryArray;

  int numLocals;
  int cluIndex;
  int numEntries = table->size();
  int *locals;

  i = 0;
  for (; index < numProcs; ++index) {
    startIndex[index] = i;
    i += numClusters[index];
  }
  totalClusters = i;

  myMinCluIndex = startIndex[myRank];

  for (index = 0; index < numLocalVertices; ++index) {
#ifdef DEBUG_COARSENER
    assert(matchVector[index] != -1);
#endif

    if (matchVector[index] < NON_LOCAL_MATCH)
      matchVector[index] += myMinCluIndex;
  }

  // ###
  // now set the non-local matches' cluster indices
  // ###

  entryArray = table->get_entries();

  for (index = 0; index < numEntries; ++index) {
    entry_ = entryArray[index];

#ifdef DEBUG_COARSENER
    assert(entry_);
#endif

    cluIndex = entry_->cluster_index();

#ifdef DEBUG_COARSENER
    assert(cluIndex >= 0);
#endif

    if (cluIndex != MATCHED_LOCALLY) {
      numLocals = entry_->number_local();
      locals = entry_->local_vertices_array();

#ifdef DEBUG_COARSENER
      assert(locals);
#endif
      cluIndex += startIndex[entry_->non_local_process()];

      // indicate - there is no reason now why a vertex may not
      // match non-locally with a vertex that was matched non-locally
      // with a vertex from another processor during a high-to-low
      // communication for example...otherwise the cluster weight
      // information is redundant?

      for (i = 0; i < numLocals; ++i)
        matchVector[locals[i] - minVertexIndex] = cluIndex;
    }
  }

#ifdef DEBUG_COARSENER
  for (index = 0; index < numLocalVertices; ++index)
    if (matchVector[index] < 0 || matchVector[index] >= totalVertices) {
      std::cout << "matchVector[" << index << "]  = " << matchVector[index] << std::endl;
      assert(0);
    }
#endif
}

int ParaFCCoarsener::accept(int locVertex, int nonLocCluWt, int highToLow,
                            int maxWt) {
  int locVertexIndex = locVertex - minVertexIndex;
  int matchValue = matchVector[locVertexIndex];

  if (matchValue < NON_LOCAL_MATCH) {
    if (clusterWeights[matchValue] + nonLocCluWt >= maxWt)
      return 0;
    else {
      // ###
      // the non-local requesting vertices will have the same
      // cluster index as locVertex
      // ###

      clusterWeights[matchValue] += nonLocCluWt;
      return 1;
    }
  } else {

    int nonLocReq = matchValue - NON_LOCAL_MATCH;
    int proc = nonLocReq / (totalVertices / numProcs);

    if ((highToLow && proc < myRank) || (!highToLow && proc > myRank))
      return 0;

    // ###
    // here think about allowing more than one cross-processor match
    // ###

    if (table->cluster_index(nonLocReq) != -1)
      return 0;

    int cluWt = vWeight[locVertexIndex] + nonLocCluWt;

    if (cluWt >= maxWt)
      return 0;

    else {
      // ###
      // the non-local requesting vertices will have the
      // same cluster index as locVertex
      // remove locVertex from the request table
      // ###

      table->remove_local(nonLocReq, locVertex, vWeight[locVertexIndex]);
      matchVector[locVertexIndex] = clusterIndex;
      clusterWeights.assign(clusterIndex++, cluWt);

      return 1;
    }
  }
}

void ParaFCCoarsener::printVisitOrder(int variable) const {
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
