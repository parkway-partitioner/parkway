
#ifndef _PARA_APPROX_FCCOARSENER_CPP
#define _PARA_APPROX_FCCOARSENER_CPP

// ### ParaApproxFCCoarsener.cpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 31/12/2004: Last Modified
//
// ###

#include "ParaApproxFCCoarsener.hpp"

ParaApproxFCCoarsener::ParaApproxFCCoarsener(int rank, int nProcs, int nParts,
                                             int percentile, int inc,
                                             int vertVisOrder,
                                             int matchReqOrder, int divByWt,
                                             int divByLen, ostream &out)
    : ParaApproxCoarsener(rank, nProcs, nParts, percentile, inc, out) {
  vertexVisitOrder = vertVisOrder;
  matchRequestVisitOrder = matchReqOrder;
  divByCluWt = divByWt;
  divByHedgeLen = divByLen;
  limitOnIndexDuringCoarsening = 0;

  table = NULL;
  // connTable = NULL;
}

ParaApproxFCCoarsener::~ParaApproxFCCoarsener() {
  DynaMem<MatchRequestTable>::deletePtr(table);
  // DynaMem<ConnVertTable>::deletePtr(connTable);
}

void ParaApproxFCCoarsener::dispCoarseningOptions() const {
  switch (dispOption) {
  case SILENT:

    break;

  default:

    out_stream << "|--- PARA_C:" << endl
               << "|- ApproxPFC:"
               << " r = " << reductionRatio << " min = " << minNodes
               << " vvo = ";
    printVisitOrder(vertexVisitOrder);
    out_stream << " mvo = ";
    printVisitOrder(matchRequestVisitOrder);
    out_stream << " divWt = " << divByCluWt << " divLen = " << divByHedgeLen
               << " %ile = " << startPercentile << " inc = " << increment
               << endl
               << "|" << endl;
    break;
  }
}

void ParaApproxFCCoarsener::buildAuxiliaryStructs(int numTotPins,
                                                  double aveVertDeg,
                                                  double aveHedgeSize) {
  // ###
  // build the MatchRequestTable
  // ###

  int i =
      static_cast<int>(ceil(static_cast<double>(numTotPins) / aveVertDeg));

  // table = new MatchRequestTable(Funct::setTableSize(i/numProcs));
  table = new MatchRequestTable(TableUtils::getTableSize(i / numProcs));
  // i = Shiftl(static_cast<int>(ceil(aveVertDeg*aveHedgeSize)),4);

  // ###
  // build the ConnVertTable
  // ###

  // connTable = new ConnVertTable(Funct::setTableSize(i));
}

void ParaApproxFCCoarsener::releaseMemory() {
  hEdgeWeight.setLength(0);
  hEdgeOffset.setLength(0);
  locPinList.setLength(0);

  vToHedgesOffset.setLength(0);
  vToHedgesList.setLength(0);
  allocHedges.setLength(0);
  clusterWeights.setLength(0);

  freeMemory();
}

ParaHypergraph *ParaApproxFCCoarsener::coarsen(ParaHypergraph &h,
                                               MPI_Comm comm) {
  loadHyperGraph(h, comm);

  if (totalVertices < minNodes || h.dontCoarsen()) {
    currPercentile = startPercentile;
    return NULL;
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

  int hEdge;
  int vertex;

  double metric;
  double maxMatchMetric;
  double reducedBy;

  MapToPosInt matchInfoLoc;

  if (numLocPins < totalVertices / 2)
    matchInfoLoc.createTable(numLocPins, 1);
  else
    matchInfoLoc.createTable(totalVertices, 0);
  // ConnVertData *vData;
  // ConnVertData **vDataArray;

  DynamicArray<int> neighVerts;
  DynamicArray<int> neighCluWts;
  DynamicArray<double> connectVals;
  DynamicArray<int> vertices(numLocalVertices);

  permuteVerticesArray(vertices.getArray(), numLocalVertices);

  if (dispOption > 1) {
    for (i = 0; i < numLocalVertices; ++i) {
      if (vWeight[i] > maxLocWt)
        maxLocWt = vWeight[i];
    }

    MPI_Reduce(&maxLocWt, &maxWt, 1, MPI_INT, MPI_MAX, 0, comm);

    if (myRank == 0) {
      out_stream << "[PFCC] " << maxVertexWt << " " << maxWt << " "
                 << aveVertexWt << " ";
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
            // ###
            // now try not to check the weight before checking the vertex
            // ###

            auto neighbourLoc = matchInfoLoc.getCareful(candidatV);
            // vData = connTable->getDataStruct(candidatV);

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
                              table->lookupClusterWt(nonLocV) + aveVertexWt;
                } else
                  cluWeight =
                      vWeight[vertex] +
                      clusterWeights[matchVector[candidatV - minVertexIndex]];
              } else {
                // ###
                // candidatV is not a local vertex or it is a local
                // vertex matched with a non-local one
                // ###

                candVwt = table->lookupClusterWt(candidatV);

                if (candVwt != -1)
                  cluWeight = vWeight[vertex] + candVwt + aveVertexWt;
                else
                  cluWeight = vWeight[vertex] + aveVertexWt;
              }

#ifdef DEBUG_COARSENER
              assert(numVisited >= 0);
#endif
              if (matchInfoLoc.insertKey(candidatV, numVisited)) {
                write_log(myRank, "numEntries %d",
                          matchInfoLoc.getNumEntries());
                write_log(myRank, "using hash %d",
                          matchInfoLoc.usingHashTable());
                write_log(myRank, "numSlots %d", matchInfoLoc.getNumSlots());
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
            /*
               if(vData)
               {
               if(divByHedgeLen)
               vData->incConNetWt(static_cast<double>(hEdgeWeight[hEdge])/(hEdgeLen-1));
               else
               vData->incConNetWt(static_cast<double>(hEdgeWeight[hEdge]));
               }
               */
            /*
               else
               {
            // ###
            // here compute the cluster weight
            // ###

            if(candidatV>=minVertexIndex && candidatV<maxVertexIndex)
            {
            // ###
            // candidatV is a local vertex
            // ###

            if(matchVector[candidatV-minVertexIndex] == -1)
            cluWeight = vWeight[vertex] + vWeight[candidatV-minVertexIndex];
            else
            if(matchVector[candidatV-minVertexIndex] >= NON_LOCAL_MATCH)
            {
            nonLocV = matchVector[candidatV-minVertexIndex] - NON_LOCAL_MATCH;
            cluWeight = vWeight[vertex] + table->lookupClusterWt(nonLocV) +
            aveVertexWt;
            }
            else
            cluWeight = vWeight[vertex] +
            clusterWeights[matchVector[candidatV-minVertexIndex]];
            }
            else
            {
            // ###
            // candidatV is not a local vertex or it is a local
            // vertex matched with a non-local one
            // ###

            candVwt = table -> lookupClusterWt(candidatV);

            if(candVwt != -1) cluWeight = vWeight[vertex] + candVwt +
            aveVertexWt;
            else cluWeight = vWeight[vertex] + aveVertexWt;
            }

            if(divByHedgeLen)
            vData = new
            ConnVertData(candidatV,static_cast<double>(hEdgeWeight[hEdge])/(hEdgeLen-1),cluWeight,NULL);
            else
            vData = new
            ConnVertData(candidatV,static_cast<double>(hEdgeWeight[hEdge]),cluWeight,NULL);

            connTable->addDataStruct(vData,candidatV);
            }

*/
          }
        }
      }

      // ###
      // now choose best match from adjacent vertices
      // visited above
      // ###

      // numVisited = connTable->getNumEntries();
      // vDataArray = connTable->getDataArray();

      // ###
      // pick best match
      // ###

      for (i = 0; i < numVisited; ++i) {
        /*
#  ifdef DEBUG_COARSENER
assert(vDataArray[i]);
#  endif
vData = vDataArray[i];
pairWt = vData->getPairWt();
*/

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
        /*
           if(pairWt <= maxVertexWt)
           {
           metric = vData->getConNetWt();

           if(divByCluWt)
           metric /= pairWt;

           if(metric > maxMatchMetric)
           {
           maxMatchMetric = metric;
           bestMatch = vData->getVertex();
           bestMatchWt = pairWt;

#  ifdef DEBUG_COARSENER
assert(bestMatch >= 0);
#  endif
}
}
*/
      }

      matchInfoLoc.resetSlots();
      numVisited = 0;
      // connTable->clearTable();

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
              table->addLocal(nonLocV, vertex + minVertexIndex, vWeight[vertex],
                              min(nonLocV / vPerProc, numProcs - 1));
#ifdef DEBUG_COARSENER
              assert(min(nonLocV / vPerProc, numProcs - 1) != myRank);
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

          table->addLocal(bestMatch, vertex + minVertexIndex, vWeight[vertex],
                          min(bestMatch / vPerProc, numProcs - 1));
#ifdef DEBUG_COARSENER
          assert(min(bestMatch / vPerProc, numProcs - 1) != myRank);
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
                    (numNotMatched + clusterIndex + table->getNumEntries());
        break;
      }
    }
  }

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

  table->clearTable();

  // ###
  // now construct the coarse hypergraph using the matching vector
  // ###

  return (contractHyperedges(h, comm));
}

void ParaApproxFCCoarsener::setRequestArrays(int highToLow) {
  int numRequests = table->getNumEntries();
  int nonLocVertex;
  int cluWt;

  int i;
  int procRank;

  MatchRequestEntry *entry;
  MatchRequestEntry **entryArray = table->getEntriesArray();

  for (i = 0; i < numProcs; ++i)
    sendLens[i] = 0;

  for (i = 0; i < numRequests; ++i) {
    entry = entryArray[i];

#ifdef DEBUG_COARSENER
    assert(entry);
#endif

    nonLocVertex = entry->getNonLocal();
    cluWt = entry->getClusterWt();
    procRank = entry->getNonLocProc();

#ifdef DEBUG_COARSENER
    assert(procRank < numProcs);
#endif

    if ((cluWt > 0) && ((highToLow && procRank < myRank) ||
                        (!highToLow && procRank > myRank))) {
      dataOutSets[procRank]->assign(sendLens[procRank]++, nonLocVertex);
      dataOutSets[procRank]->assign(sendLens[procRank]++, cluWt);
    }
  }
}

void ParaApproxFCCoarsener::setReplyArrays(int highToLow, int maxVWt) {
  int j;
  int i;
  int l;

  DynamicArray<int> visitOrder;

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
      visitOrder.setLength(visitOrderLen);

      for (j = 0; j < visitOrderLen; ++j)
        visitOrder[j] = j;

      // ###
      // processing match requests in random order
      // ###

      Funct::randomPermutation(visitOrder.getArray(), visitOrderLen);

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

void ParaApproxFCCoarsener::processReqReplies() {
  int i;
  int j;
  int index;

  int startOffset = 0;
  int vNonLocReq;
  int cluWt;
  int matchIndex;
  int numLocals;
  int *locals;

  MatchRequestEntry *entry;

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
        entry = table->getEntryPtr(vNonLocReq);

#ifdef DEBUG_COARSENER
        assert(entry);
#endif

        entry->setCluIndex(matchIndex);
        entry->setCluWeight(cluWt);
      } else {
        // ###
        // match not successful - match requesting
        // vertices into a cluster
        // ###

        entry = table->getEntryPtr(vNonLocReq);
        locals = entry->getLocalsArray();
        numLocals = entry->getNumLocals();
        entry->setCluIndex(MATCHED_LOCALLY);

        for (index = 0; index < numLocals; ++index)
          matchVector[locals[index] - minVertexIndex] = clusterIndex;

        clusterWeights.assign(clusterIndex++, entry->getClusterWt());
      }
    }
    startOffset += recvLens[i];
  }
}

void ParaApproxFCCoarsener::permuteVerticesArray(int *verts, int nLocVerts) {
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

void ParaApproxFCCoarsener::setClusterIndices(MPI_Comm comm) {
  DynamicArray<int> numClusters(numProcs);
  DynamicArray<int> startIndex(numProcs);

  MPI_Allgather(&clusterIndex, 1, MPI_INT, numClusters.getArray(), 1, MPI_INT,
                comm);

  int index = 0;
  int i;

  MatchRequestEntry *entry;
  MatchRequestEntry **entryArray;

  int numLocals;
  int cluIndex;
  int numEntries = table->getNumEntries();
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

  entryArray = table->getEntriesArray();

  for (index = 0; index < numEntries; ++index) {
    entry = entryArray[index];

#ifdef DEBUG_COARSENER
    assert(entry);
#endif

    cluIndex = entry->getCluIndex();

    if (cluIndex >= 0 && cluIndex != MATCHED_LOCALLY) {
      numLocals = entry->getNumLocals();
      locals = entry->getLocalsArray();

#ifdef DEBUG_COARSENER
      assert(locals);
#endif
      cluIndex += startIndex[entry->getNonLocProc()];

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
    assert(matchVector[index] >= 0 && matchVector[index] < totalVertices);
#endif
}

int ParaApproxFCCoarsener::accept(int locVertex, int nonLocCluWt, int highToLow,
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

    if (table->lookupCluIndex(nonLocReq) != -1)
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

      table->removeLocal(nonLocReq, locVertex, vWeight[locVertexIndex]);
      matchVector[locVertexIndex] = clusterIndex;
      clusterWeights.assign(clusterIndex++, cluWt);

      return 1;
    }
  }
}

void ParaApproxFCCoarsener::printVisitOrder(int variable) const {
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
