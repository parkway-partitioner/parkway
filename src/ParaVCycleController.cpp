
#ifndef _VCYCLE_CONTROLLER_CPP
#define _VCYCLE_CONTROLLER_CPP

// ### ParaVCycleController.cpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 4/1/2005: Last Modified
//
// ###

#include "ParaVCycleController.hpp"

ParaVCycleController::ParaVCycleController(ParaRestrCoarsener &rc,
                                           ParaCoarsener &c, ParaRefiner &r,
                                           SeqController &ref, int rank, int nP,
                                           int percentile, int inc,
                                           int approxRef, int limit,
                                           double limitAsPercent, ostream &out)
    : ParaController(c, r, ref, rank, nP, percentile, inc, approxRef, out),
      restrCoarsener(rc) {
  if (limit)
    limitOnCycles = limit;
  else
    limitOnCycles = LARGE_CONSTANT;

  if (limitAsPercent)
    limAsPercentOfCut = limitAsPercent;
  else
    limitAsPercent = 0;

  minInterVertIndex = 0;

  mapToInterVerts.setLength(0);
  mapToOrigVerts.setLength(0);
}

ParaVCycleController::~ParaVCycleController() {}

void ParaVCycleController::dispParaControllerOptions() const {
  switch (dispOption) {
  case SILENT:
    break;

  default:

    out_stream << "|--- PARA_CONTR (# parts = " << numTotalParts
               << "): " << endl
               << "|- VCYCLE:"
               << " pRuns = " << numParaRuns << " kT = " << keepPartitionsWithin
               << " rKT = " << reductionInKeepThreshold
               << " appRef = " << approxRefine
               << " wTF = " << writePartitionToFile;
    printVCycleType();
    out_stream << " lim = " << limitOnCycles << " %min = " << limAsPercentOfCut
               << " start %le = " << startPercentile
               << " %le inc = " << percentileIncrement << endl
               << "|" << endl;
    break;
  }
}

void ParaVCycleController::setWeightConstraints(MPI_Comm comm) {
#ifdef DEBUG_CONTROLLER
  assert(hgraph);
  assert(numTotalParts != -1);
  assert(balConstraint > 0.0 && balConstraint < 1.0);
#endif

  int locGraphWt;
  int totGraphWt;
  int maxVertWt;

  double avePartWt;

  locGraphWt = hgraph->getLocalVertexWt();

  MPI_Allreduce(&locGraphWt, &totGraphWt, 1, MPI_INT, MPI_SUM, comm);

  avePartWt = static_cast<double>(totGraphWt) / numTotalParts;
  maxPartWt = static_cast<int>(floor(avePartWt + avePartWt * balConstraint));
  maxVertWt = static_cast<int>(floor(avePartWt * balConstraint));

  coarsener.setMaxVertexWt(maxVertWt);
  restrCoarsener.setMaxVertexWt(maxVertWt);

  coarsener.setTotGraphWt(totGraphWt);
  restrCoarsener.setTotGraphWt(totGraphWt);

  seqController.setMaxVertexWt(maxVertWt);
}

void ParaVCycleController::recordVCyclePartition(const ParaHypergraph &h,
                                                 int numIteration) {
#ifdef DEBUG_CONTROLLER
  assert(numIteration >= 0);
#endif

  int i;

  int minVertexIndex;
  int numLocalVertices;

  int *array;
  int *pVector;

  IntArray *bestPartVector;

  pVector = h.getPartVectorArray();
  minVertexIndex = h.getMinVertexIndex();
  numLocalVertices = h.getNumLocalVertices();

  if (numIteration == 0) {
    bestPartVector = new IntArray(numLocalVertices);

    array = bestPartVector->getArray();

    for (i = 0; i < numLocalVertices; ++i)
      array[i] = pVector[i];

    minLocCurrVertId.push(minVertexIndex);
    numLocCurrVerts.push(numLocalVertices);
    bestVCyclePartition.push(bestPartVector);
  } else {
    bestPartVector = bestVCyclePartition.getTopElem();

#ifdef DEBUG_CONTROLLER
    assert(bestPartVector);
#endif

    bestPartVector->setLength(numLocalVertices);
    array = bestPartVector->getArray();

    for (i = 0; i < numLocalVertices; ++i)
      array[i] = pVector[i];

    minLocCurrVertId.pop();
    numLocCurrVerts.pop();
    minLocCurrVertId.push(minVertexIndex);
    numLocCurrVerts.push(numLocalVertices);
  }

#ifdef DEBUG_CONTROLLER
  assert(minLocCurrVertId.getNumElem() > 0);
  assert(numLocCurrVerts.getNumElem() > 0);
  assert(h.getNumPartitions() == 1);
  assert(bestVCyclePartition.getNumElem() > 0);
#endif
}

void ParaVCycleController::gatherInVCyclePartition(ParaHypergraph &h, int cut,
                                                   MPI_Comm comm) {
#ifdef DEBUG_CONTROLLER
  assert(minLocCurrVertId.getNumElem() > 0);
  assert(numLocCurrVerts.getNumElem() > 0);
  assert(h.getNumPartitions() == 1);
  assert(bestVCyclePartition.getNumElem() > 0);
#endif

  int i;
  int j;
  int ij;

  int minStoredVertexIndex;
  int maxStoredVertexIndex;
  int numStoredLocVertices;
  int numLocalVertices;
  int origV;
  int arrayLen;
  int totToSend;
  int totToRecv;

#ifdef DEBUG_CONTROLLER
  int numTotVertices = h.getNumTotalVertices();
#endif

  int *array;
  int *dataOutArray;
  int *myPartVector;
  int *myVtoOrigV;

  myPartVector = h.getPartVectorArray();
  numLocalVertices = h.getNumLocalVertices();
  myVtoOrigV = h.getToOrigVArray();

  minStoredVertexIndex = minLocCurrVertId.pop();
  numStoredLocVertices = numLocCurrVerts.pop();
  maxStoredVertexIndex = minStoredVertexIndex + numStoredLocVertices;

  IntArray *bestPartVector = bestVCyclePartition.pop();

#ifdef DEBUG_CONTROLLER
  assert(bestPartVector);
#endif

  array = bestPartVector->getArray();

  DynamicArray<int> minStoredIDs(numProcs);
  DynamicArray<int> localSendArrayVerts;
  DynamicArray<int> procs(numProcs);

  MPI_Allgather(&minStoredVertexIndex, 1, MPI_INT, minStoredIDs.getArray(), 1,
                MPI_INT, comm);

  for (i = 0; i < numProcs; ++i)
    procs[i] = i;

  CompleteBinaryTree<int> storedVtoProc(procs.getArray(),
                                        minStoredIDs.getArray(), numProcs);

  for (i = 0; i < numProcs; ++i)
    sendLens[i] = 0;

  for (i = 0; i < numLocalVertices; ++i) {
    origV = myVtoOrigV[i];
#ifdef DEBUG_CONTROLLER
    assert(origV >= 0 && origV < numTotVertices);
#endif

    if (origV >= minStoredVertexIndex && origV < maxStoredVertexIndex) {
      myPartVector[i] = array[origV - minStoredVertexIndex];
    } else {
      j = storedVtoProc.getRootVal(origV);

      dataOutSets[j]->assign(sendLens[j]++, origV);
      dataOutSets[j]->assign(sendLens[j]++, i);
    }
  }

  ij = 0;

  for (i = 0; i < numProcs; ++i) {
    sendLens[i] = Shiftr(sendLens[i], 1);
    sendDispls[i] = ij;
    ij += sendLens[i];
  }

  sendArray.setLength(ij);
  localSendArrayVerts.setLength(ij);
  totToSend = ij;

  j = 0;

  for (i = 0; i < numProcs; ++i) {
    dataOutArray = dataOutSets[i]->getArray();
    arrayLen = Shiftl(sendLens[i], 1);

    for (ij = 0; ij < arrayLen;) {
      sendArray[j] = dataOutArray[ij++];
      localSendArrayVerts[j++] = dataOutArray[ij++];
    }
  }
#ifdef DEBUG_CONTROLLER
  assert(j == totToSend);
#endif

  MPI_Alltoall(sendLens.getArray(), 1, MPI_INT, recvLens.getArray(), 1, MPI_INT,
               comm);

  j = 0;
  for (i = 0; i < numProcs; ++i) {
    recvDispls[i] = j;
    j += recvLens[i];
  }

  receiveArray.setLength(j);
  totToRecv = j;

  MPI_Alltoallv(sendArray.getArray(), sendLens.getArray(),
                sendDispls.getArray(), MPI_INT, receiveArray.getArray(),
                recvLens.getArray(), recvDispls.getArray(), MPI_INT, comm);

  // ###
  // now have received all requests and sent out our requests
  // the reply communication will have the dual dimensions
  // ###

  sendArray.setLength(totToRecv);

  for (i = 0; i < totToRecv; ++i) {
#ifdef DEBUG_CONTROLLER
    assert(receiveArray[i] >= minStoredVertexIndex &&
           receiveArray[i] < maxStoredVertexIndex);
#endif

    sendArray[i] = array[receiveArray[i] - minStoredVertexIndex];
  }

  receiveArray.setLength(totToSend);

  MPI_Alltoallv(sendArray.getArray(), recvLens.getArray(),
                recvDispls.getArray(), MPI_INT, receiveArray.getArray(),
                sendLens.getArray(), sendDispls.getArray(), MPI_INT, comm);

  for (i = 0; i < totToSend; ++i) {
    myPartVector[localSendArrayVerts[i]] = receiveArray[i];
  }

  h.setCut(0, cut);

#ifdef DEBUG_CONTROLLER
  h.checkPartitions(numTotalParts, maxPartWt, comm);
#endif

  DynaMem<IntArray>::deletePtr(bestPartVector);
}

void ParaVCycleController::projectVCyclePartition(ParaHypergraph &cG,
                                                  ParaHypergraph &fG,
                                                  MPI_Comm comm) {
  // function spec:
  // if mapToInterVerts is empty (i.e. the map to the
  // then create a default map (i.e. map[0] = 0 etc.)
  // do a normal projection
  // else
  // do a special projection where the partition id projected
  // via the mapToInterVerts
  // create a default map (i.e. map[0] = 0 etc.)
  // do a normal projection

  int i;

  int numLocalFineVertices = fG.getNumLocalVertices();

  if (mapToInterVerts.getLength() != 0) {
    fG.setNumberPartitions(1);

    int numLocalCoarseVertices = cG.getNumLocalVertices();
    int numTotalCoarseVertices = cG.getNumTotalVertices();
    int coarseCut = cG.getCut(0);

    fG.setCut(0, coarseCut);

#ifdef DEBUG_CONTROLLER
    assert(numLocalCoarseVertices == mapToInterVerts.getLength());
    if (myRank != numProcs - 1)
      assert(numLocalCoarseVertices == numTotalCoarseVertices / numProcs);
#endif

    int *coarsePartVector = cG.getPartVectorArray();
    int *finePartVector = fG.getPartVectorArray();
    int *fineMatVector = fG.getMatchVectorArray();
    int *array;

    int minCoarseVertexId = cG.getMinVertexIndex();
#ifdef DEBUG_CONTROLLER
    int maxCoarseVertexId = minCoarseVertexId + numLocalCoarseVertices;
#endif
    int numRequestingLocalVerts;
    int cVertPerProc = numTotalCoarseVertices / numProcs;
    int cVertex;
    int vPart;
    int totToRecv;
    int totToSend;
    int sendLength;

    int j;
    int ij;

    DynamicArray<int> interGraphPVector(numLocalCoarseVertices);
    DynamicArray<int> requestingLocalVerts;

    for (i = 0; i < numProcs; ++i)
      sendLens[i] = 0;

    for (i = 0; i < numLocalCoarseVertices; ++i) {
      cVertex = mapToInterVerts[i];
      vPart = coarsePartVector[i];

#ifdef DEBUG_CONTROLLER
      assert(cVertex >= 0 && cVertex < numTotalCoarseVertices);
      assert(vPart >= 0 && vPart < numTotalParts);
#endif

      ij = min(cVertex / cVertPerProc, numProcs - 1);

      if (ij == myRank) {
        interGraphPVector[cVertex - minCoarseVertexId] = vPart;
      } else {
        dataOutSets[ij]->assign(sendLens[ij]++, cVertex);
        dataOutSets[ij]->assign(sendLens[ij]++, vPart);
      }
    }

    ij = 0;

    for (i = 0; i < numProcs; ++i) {
      sendDispls[i] = ij;
      ij += sendLens[i];
    }

    sendArray.setLength(ij);
    totToSend = ij;
    ij = 0;

    for (i = 0; i < numProcs; ++i) {
      j = 0;
      sendLength = sendLens[i];
      array = dataOutSets[i]->getArray();

      while (j < sendLength) {
        sendArray[ij++] = array[j++];
      }
    }

#ifdef DEBUG_CONTROLLER
    assert(ij == totToSend);
#endif

    // ###
    // get dimension and carry
    // out the communication
    // ###

    MPI_Alltoall(sendLens.getArray(), 1, MPI_INT, recvLens.getArray(), 1,
                 MPI_INT, comm);

    ij = 0;
    for (i = 0; i < numProcs; ++i) {
      recvDispls[i] = ij;
      ij += recvLens[i];
    }

    receiveArray.setLength(ij);
    totToRecv = ij;

    MPI_Alltoallv(sendArray.getArray(), sendLens.getArray(),
                  sendDispls.getArray(), MPI_INT, receiveArray.getArray(),
                  recvLens.getArray(), recvDispls.getArray(), MPI_INT, comm);

    // ###
    // now finish the initialisation
    // of the corse hypergraph partitions
    // ###

    for (i = 0; i < totToRecv;) {
      cVertex = receiveArray[i++];
      vPart = receiveArray[i++];

#ifdef DEBUG_CONTROLLER
      assert(cVertex >= minCoarseVertexId && cVertex < maxCoarseVertexId);
      assert(vPart >= 0 && vPart < numTotalParts);
#endif
      interGraphPVector[cVertex - minCoarseVertexId] = vPart;
    }

#ifdef DEBUG_CONTROLLER
    for (i = 0; i < numLocalCoarseVertices; ++i)
      assert(interGraphPVector[i] >= 0 && interGraphPVector[i] < numTotalParts);
#endif

    // ###
    // now have the coarse hypergraph partition initialised
    // ###

    for (i = 0; i < numProcs; ++i)
      sendLens[i] = 0;

    // ###
    // initialise the local partition structures
    // ###

    for (i = 0; i < numLocalFineVertices; ++i) {
#ifdef DEBUG_CONTROLLER
      assert(fineMatVector[i] >= 0 &&
             fineMatVector[i] < numTotalCoarseVertices);
#endif
      cVertex = fineMatVector[i];
      ij = min(cVertex / cVertPerProc, numProcs - 1);

      if (ij == myRank) {
        finePartVector[i] = interGraphPVector[cVertex - minCoarseVertexId];
      } else {
        dataOutSets[ij]->assign(sendLens[ij]++, i);
        dataOutSets[ij]->assign(sendLens[ij]++, cVertex);
      }
    }

    // ###
    // prepare to send requests for partition vector values
    // ###

    ij = 0;

    for (i = 0; i < numProcs; ++i) {
      sendDispls[i] = ij;
      ij += (Shiftr(sendLens[i], 1));
    }

    sendArray.setLength(ij);
    requestingLocalVerts.setLength(ij);

    numRequestingLocalVerts = ij;
    totToSend = ij;
    ij = 0;

    for (i = 0; i < numProcs; ++i) {
      j = 0;
      sendLength = sendLens[i];
      array = dataOutSets[i]->getArray();

      while (j < sendLength) {
        requestingLocalVerts[ij] = array[j++];
        sendArray[ij++] = array[j++];
      }

      sendLens[i] = Shiftr(sendLength, 1);
    }

#ifdef DEBUG_CONTROLLER
    assert(ij == totToSend);
#endif

    // ###
    // get dimension and carry
    // out the communication
    // ###

    MPI_Alltoall(sendLens.getArray(), 1, MPI_INT, recvLens.getArray(), 1,
                 MPI_INT, comm);

    ij = 0;
    for (i = 0; i < numProcs; ++i) {
      recvDispls[i] = ij;
      ij += recvLens[i];
    }

    receiveArray.setLength(ij);
    totToRecv = ij;

    MPI_Alltoallv(sendArray.getArray(), sendLens.getArray(),
                  sendDispls.getArray(), MPI_INT, receiveArray.getArray(),
                  recvLens.getArray(), recvDispls.getArray(), MPI_INT, comm);

    // ###
    // process the requests for
    // local vertex partitions
    // dimensions of cummunication
    // are reversed!
    // ###

    totToSend = totToRecv;
    sendArray.setLength(totToSend);

    for (i = 0; i < totToRecv; ++i) {
#ifdef DEBUG_CONTROLLER
      assert(receiveArray[i] >= minCoarseVertexId &&
             receiveArray[i] < maxCoarseVertexId);
#endif
      cVertex = receiveArray[i] - minCoarseVertexId;
      sendArray[i] = interGraphPVector[cVertex];
    }

    totToRecv = numRequestingLocalVerts;
    receiveArray.setLength(totToRecv);

    // ###
    // now do the dual communication
    // ###

    MPI_Alltoallv(sendArray.getArray(), recvLens.getArray(),
                  recvDispls.getArray(), MPI_INT, receiveArray.getArray(),
                  sendLens.getArray(), sendDispls.getArray(), MPI_INT, comm);

    // ###
    // finish off initialising the
    // partition vector
    // ###

    ij = 0;

    for (i = 0; i < numRequestingLocalVerts; ++i) {
      cVertex = requestingLocalVerts[i];
#ifdef DEBUG_CONTROLLER
      assert(cVertex >= 0 && cVertex < numLocalFineVertices);
      assert(receiveArray[ij] >= 0 && receiveArray[ij] < numTotalParts);
#endif
      finePartVector[cVertex] = receiveArray[ij++];
    }

#ifdef DEBUG_CONTROLLER
    for (i = 0; i < numLocalFineVertices; ++i)
      assert(finePartVector[i] >= 0 && finePartVector[i] < numTotalParts);
    int calcCut = fG.calcCutsize(numTotalParts, 0, comm);
    assert(coarseCut == calcCut);
#endif
  } else {
    fG.projectPartitions(cG, comm);
  }

  mapToInterVerts.setLength(numLocalFineVertices);
  minInterVertIndex = fG.getMinVertexIndex();

  for (i = 0; i < numLocalFineVertices; ++i)
    mapToInterVerts[i] = minInterVertIndex + i;
}

void ParaVCycleController::shuffleVCycleVertsByPartition(ParaHypergraph &h,
                                                         MPI_Comm comm) {
  // function spec:
  // as the vertices are shuffled modify the mapToInterVerts
  // check that the sum of the lengths of the mapToInterVerts
  // across the processors is equal to the number of totalVertices
  // of h
  //

  int i;
  int j;
  int ij;

  int minLocVertIdBefShuff = h.getMinVertexIndex();

#ifdef DEBUG_CONTROLLER
  int numLocVertBefShuff = h.getNumLocalVertices();
  int maxLocVertIdBefShuff = minLocVertIdBefShuff + numLocVertBefShuff;
#endif

  int numTotVertices = h.getNumTotalVertices();
  int numVBefPerProc = numTotVertices / numProcs;

#ifdef DEBUG_CONTROLLER
  assert(minInterVertIndex == minLocVertIdBefShuff);
  assert(numLocVertBefShuff == mapToInterVerts.getLength());
#endif

  h.shuffleVerticesByPartition(numTotalParts, comm);

  int numLocVertAftShuff = h.getNumLocalVertices();
  int minLocVertIdAftShuff = h.getMinVertexIndex();
  int totToSend;
  int totToRecv;
  int sendLength;
  int numRequestingLocalVerts;

  int *vToOrigV = h.getToOrigVArray();
  int *array;

  DynamicArray<int> *newToInter = new DynamicArray<int>(numLocVertAftShuff);
  DynamicArray<int> requestingLocalVerts;

  for (i = 0; i < numProcs; ++i)
    sendLens[i] = 0;

  array = newToInter->getArray();

  for (i = 0; i < numLocVertAftShuff; ++i) {
    j = vToOrigV[i];
    ij = min(j / numVBefPerProc, numProcs - 1);

    if (ij == myRank) {
      array[i] = mapToInterVerts[j - minInterVertIndex];
    } else {
      dataOutSets[ij]->assign(sendLens[ij]++, i);
      dataOutSets[ij]->assign(sendLens[ij]++, j);
    }
  }

  ij = 0;
  for (i = 0; i < numProcs; ++i) {
    sendDispls[i] = ij;
    ij += (Shiftr(sendLens[i], 1));
  }

  sendArray.setLength(ij);
  requestingLocalVerts.setLength(ij);
  numRequestingLocalVerts = ij;
  totToSend = ij;
  ij = 0;

  for (i = 0; i < numProcs; ++i) {
    j = 0;
    sendLength = sendLens[i];
    array = dataOutSets[i]->getArray();

    while (j < sendLength) {
      requestingLocalVerts[ij] = array[j++];
      sendArray[ij++] = array[j++];
    }

    sendLens[i] = Shiftr(sendLens[i], 1);
  }

#ifdef DEBUG_CONTROLLER
  assert(ij == totToSend);
#endif

  // ###
  // get dimension and carry
  // out the communication
  // ###

  MPI_Alltoall(sendLens.getArray(), 1, MPI_INT, recvLens.getArray(), 1, MPI_INT,
               comm);

  ij = 0;
  for (i = 0; i < numProcs; ++i) {
    recvDispls[i] = ij;
    ij += recvLens[i];
  }

  receiveArray.setLength(ij);
  totToRecv = ij;

  MPI_Alltoallv(sendArray.getArray(), sendLens.getArray(),
                sendDispls.getArray(), MPI_INT, receiveArray.getArray(),
                recvLens.getArray(), recvDispls.getArray(), MPI_INT, comm);

  // ###
  // process the requests for
  // local mapToInterVerts
  // dimensions of cummunication
  // are reversed!
  // ###

  totToSend = totToRecv;
  sendArray.setLength(totToSend);

  for (i = 0; i < totToRecv; ++i) {
#ifdef DEBUG_CONTROLLER
    assert(receiveArray[i] >= minLocVertIdBefShuff &&
           receiveArray[i] < maxLocVertIdBefShuff);
#endif
    j = receiveArray[i] - minLocVertIdBefShuff;
    sendArray[i] = mapToInterVerts[j];
  }

  totToRecv = numRequestingLocalVerts;
  receiveArray.setLength(totToRecv);

  // ###
  // now do the dual communication
  // ###

  MPI_Alltoallv(sendArray.getArray(), recvLens.getArray(),
                recvDispls.getArray(), MPI_INT, receiveArray.getArray(),
                sendLens.getArray(), sendDispls.getArray(), MPI_INT, comm);

  // ###
  // finish off initialising the
  // new mapToInterVerts vector
  // ###

  array = newToInter->getArray();

  for (i = 0; i < numRequestingLocalVerts; ++i) {
    j = requestingLocalVerts[i];
#ifdef DEBUG_CONTROLLER
    assert(j >= 0 && j < numLocVertAftShuff);
    assert(receiveArray[i] >= 0 && receiveArray[i] < numTotVertices);
#endif
    array[j] = receiveArray[i];
  }

  minInterVertIndex = minLocVertIdAftShuff;
  mapToInterVerts.setArray(array, numLocVertAftShuff);
}

/*
void ParaVCycleController::randomVCycleVertShuffle(ParaHypergraph &h,
ParaHypergraph &fineH, MPI_Comm comm)
{

  int i;
  int j;
  int ij;

  //int minLocVertIdBefShuff = h.getMinVertexIndex();

  //#  ifdef DEBUG_CONTROLLER
  //int numLocVertBefShuff = h.getNumLocalVertices();
  //int maxLocVertIdBefShuff =  minLocVertIdBefShuff+ numLocVertBefShuff;
  //#  endif

  int numTotVertices = h.getNumTotalVertices();
  //int numVBefPerProc = numTotVertices / numProcs;
  int numVPerProc = numTotVertices / numProcs;
  int minLocVertexIndex = h.getMinVertexIndex();
  int numLocVertices = h.getNumLocalVertices();

#  ifdef DEBUG_CONTROLLER
  int maxLocVertexIndex = minLocVertexIndex+numLocVertices;
#  endif

#  ifdef DEBUG_CONTROLLER
  assert(numLocVertices == mapToInterVerts.getLength());
#  endif

  //h.shuffleVerticesByPartition(numTotalParts,comm);
  h.randomVertexShuffle(fineH, comm);

  //int numLocVertAftShuff = h.getNumLocalVertices();
  //int minLocVertIdAftShuff = h.getMinVertexIndex();
  int totToSend;
  int totToRecv;
  int sendLength;
  int numRequestingLocalVerts;

  int *vToOrigV = h.getToOrigVArray();
  int *array;

  DynamicArray<int>* newToInter = new DynamicArray<int>(numLocVertices);
  DynamicArray<int> requestingLocalVerts;

  for (i=0;i<numProcs;++i)
    sendLens[i] = 0;

  array = newToInter->getArray();

  for (i=0;i<numLocVertices;++i)
    {
      j = vToOrigV[i];
      ij = min(j / numVPerProc, numProcs-1);

      if(ij == myRank)
        {
          array[i] = mapToInterVerts[j-minInterVertIndex];
        }
      else
        {
          dataOutSets[ij]->assign(sendLens[ij]++, i);
          dataOutSets[ij]->assign(sendLens[ij]++, j);
        }
    }

  ij = 0;
  for (i=0;i<numProcs;++i)
    {
      sendDispls[i] = ij;
      ij += (Shiftr(sendLens[i],1));
    }

  sendArray.setLength(ij);
  requestingLocalVerts.setLength(ij);
  numRequestingLocalVerts = ij;
  totToSend = ij;
  ij = 0;

  for (i=0;i<numProcs;++i)
    {
      j=0;
      sendLength = sendLens[i];
      array = dataOutSets[i]->getArray();

      while(j < sendLength)
        {
          requestingLocalVerts[ij] = array[j++];
          sendArray[ij++] = array[j++];
        }

      sendLens[i] = Shiftr(sendLens[i],1);
    }

#  ifdef DEBUG_CONTROLLER
  assert(ij == totToSend);
#  endif

  // ###
  // get dimension and carry
  // out the communication
  // ###

  MPI_Alltoall(sendLens.getArray(), 1, MPI_INT, recvLens.getArray(), 1, MPI_INT,
comm);

  ij = 0;
  for (i=0;i<numProcs;++i)
    {
      recvDispls[i] = ij;
      ij += recvLens[i];
    }

  receiveArray.setLength(ij);
  totToRecv = ij;

  MPI_Alltoallv(sendArray.getArray(), sendLens.getArray(),
sendDispls.getArray(), MPI_INT, receiveArray.getArray(), recvLens.getArray(),
recvDispls.getArray(), MPI_INT, comm);

  // ###
  // process the requests for
  // local mapToInterVerts
  // dimensions of cummunication
  // are reversed!
  // ###

  totToSend = totToRecv;
  sendArray.setLength(totToSend);

  for (i=0;i<totToRecv;++i)
    {
#  ifdef DEBUG_CONTROLLER
      assert(receiveArray[i] >= minLocVertexIndex && receiveArray[i] <
maxLocVertexIndex);
#  endif
      j = receiveArray[i]-minLocVertexIndex;//minLocVertIdBefShuff;
      sendArray[i] = mapToInterVerts[j];
    }

  totToRecv = numRequestingLocalVerts;
  receiveArray.setLength(totToRecv);

  // ###
  // now do the dual communication
  // ###

  MPI_Alltoallv(sendArray.getArray(), recvLens.getArray(),
recvDispls.getArray(), MPI_INT, receiveArray.getArray(), sendLens.getArray(),
sendDispls.getArray(), MPI_INT, comm);

  // ###
  // finish off initialising the
  // new mapToInterVerts vector
  // ###

  array = newToInter->getArray();

  for (i=0;i<numRequestingLocalVerts;++i)
    {
      j = requestingLocalVerts[i];
#  ifdef DEBUG_CONTROLLER
      assert(j >= 0 && j < numLocVertices);
      assert(receiveArray[i] >= 0 && receiveArray[i] < numTotVertices);
#  endif
      array[j] = receiveArray[i];
    }

  minInterVertIndex = minLocVertexIndex;
  mapToInterVerts.setArray(array,numLocVertices);
}
*/

void ParaVCycleController::shiftVCycleVertsToBalance(ParaHypergraph &h,
                                                     MPI_Comm comm) {
  // function spec:
  // as the vertices are shuffled modify the mapToInterVerts
  // check that the sum of the lengths of the mapToInterVerts
  // across the processors is equal to the number of totalVertices
  // of h
  //

  int i;
  int j;

  int numLocalVertices = h.getNumLocalVertices();
  int minVertexIndex = h.getMinVertexIndex();
  int maxVertexIndex = minVertexIndex + numLocalVertices;
  int numTotVertices = h.getNumTotalVertices();
  int vPerProc = numTotVertices / numProcs;
  int numMyNewVertices;

  if (myRank != numProcs - 1)
    numMyNewVertices = vPerProc;
  else
    numMyNewVertices = vPerProc + Mod(numTotVertices, numProcs);

  DynamicArray<int> minNewIndex(numProcs);
  DynamicArray<int> maxNewIndex(numProcs);

#ifdef DEBUG_CONTROLLER
  assert(minVertexIndex == minInterVertIndex);
#endif

  h.shiftVerticesToBalance(comm);

  // ###
  // mapToInterVerts = sendArray
  // newToInter = receiveArray
  // ###

  DynamicArray<int> *newToInter = new DynamicArray<int>(numMyNewVertices);

  for (i = 0; i < numProcs; ++i) {
    if (i == 0) {
      minNewIndex[i] = 0;
      maxNewIndex[i] = vPerProc;
    } else {
      minNewIndex[i] = maxNewIndex[i - 1];
      if (i == numProcs - 1)
        maxNewIndex[i] = numTotVertices;
      else
        maxNewIndex[i] = minNewIndex[i] + vPerProc;
    }
  }

  for (i = 0; i < numProcs; ++i) {
    if (i == 0)
      sendDispls[i] = 0;
    else
      sendDispls[i] = sendDispls[i - 1] + sendLens[i - 1];

    sendLens[i] =
        max(numLocalVertices - (max(maxVertexIndex - maxNewIndex[i], 0) +
                                max(minNewIndex[i] - minVertexIndex, 0)),
            0);
  }

  MPI_Alltoall(sendLens.getArray(), 1, MPI_INT, recvLens.getArray(), 1, MPI_INT,
               comm);

  j = 0;
  for (i = 0; i < numProcs; ++i) {
    recvDispls[i] = j;
    j += recvLens[i];
  }
#ifdef DEBUG_CONTROLLER
  assert(j == numMyNewVertices);
#endif

  MPI_Alltoallv(mapToInterVerts.getArray(), sendLens.getArray(),
                sendDispls.getArray(), MPI_INT, newToInter->getArray(),
                recvLens.getArray(), recvDispls.getArray(), MPI_INT, comm);

  minInterVertIndex = minNewIndex[myRank];
  mapToInterVerts.setArray(newToInter->getArray(), numMyNewVertices);
}

void ParaVCycleController::updateMapToOrigVerts(MPI_Comm comm) {
#ifdef DEBUG_CONTROLLER
  int numLocalVertices = hgraph->getNumLocalVertices();
#endif

  int minLocVertIndex = hgraph->getMinVertexIndex();
  int numTotalVertices = hgraph->getNumTotalVertices();
  int vertPerProc = numTotalVertices / numProcs;
  int totToRecv;
  int totToSend;
  int sendLength;
  int vertex;

  int i;
  int j;
  int ij;

#ifdef DEBUG_CONTROLLER
  assert(mapToOrigVerts.getLength() == numOrigLocVerts);
  assert(mapToInterVerts.getLength() == numOrigLocVerts);
  if (myRank != numProcs - 1)
    assert(numLocalVertices == vertPerProc);
#endif

  DynamicArray<int> *newMapToOrig = new DynamicArray<int>(numOrigLocVerts);
  DynamicArray<int> copyOfSendArray;

  int *auxArray = newMapToOrig->getArray();
  int *array;

  for (i = 0; i < numProcs; ++i)
    sendLens[i] = 0;

  for (i = 0; i < numOrigLocVerts; i++) {
    vertex = mapToInterVerts[i];

#ifdef DEBUG_CONTROLLER
    assert(vertex >= 0 && vertex < numTotalVertices);
#endif

    ij = min(vertex / vertPerProc, numProcs - 1);

    if (ij == myRank) {
      auxArray[i] = mapToOrigVerts[vertex - minLocVertIndex];
    } else {
      dataOutSets[ij]->assign(sendLens[ij]++, i);
      dataOutSets[ij]->assign(sendLens[ij]++, vertex);
    }
  }

  ij = 0;

  for (i = 0; i < numProcs; ++i) {
    sendDispls[i] = ij;
    ij += Shiftr(sendLens[i], 1);
  }

  sendArray.setLength(ij);
  copyOfSendArray.setLength(ij);
  totToSend = ij;
  ij = 0;

  for (i = 0; i < numProcs; ++i) {
    j = 0;
    sendLength = sendLens[i];
    array = dataOutSets[i]->getArray();

    while (j < sendLength) {
      copyOfSendArray[ij] = array[j++];
      sendArray[ij++] = array[j++];
    }

    sendLens[i] = Shiftr(sendLens[i], 1);
  }

#ifdef DEBUG_CONTROLLER
  assert(ij == totToSend);
#endif

  // ###
  // get dimension and carry
  // out the communication
  // ###

  MPI_Alltoall(sendLens.getArray(), 1, MPI_INT, recvLens.getArray(), 1, MPI_INT,
               comm);

  ij = 0;
  for (i = 0; i < numProcs; ++i) {
    recvDispls[i] = ij;
    ij += recvLens[i];
  }

  receiveArray.setLength(ij);
  totToRecv = ij;

  MPI_Alltoallv(sendArray.getArray(), sendLens.getArray(),
                sendDispls.getArray(), MPI_INT, receiveArray.getArray(),
                recvLens.getArray(), recvDispls.getArray(), MPI_INT, comm);

  // ###
  // now initialise the bestPartition array
  // using the data in receiveArray
  // ###

  sendArray.setLength(totToRecv);

  for (i = 0; i < totToRecv; ++i) {
#ifdef DEBUG_CONTROLLER
    assert(receiveArray[i] >= minLocVertIndex &&
           receiveArray[i] < minLocVertIndex + numOrigLocVerts);
#endif
    sendArray[i] = mapToOrigVerts[receiveArray[i] - minLocVertIndex];
  }

  receiveArray.setLength(totToSend);

  // ###
  // now do 'dual' communication
  // ###

  MPI_Alltoallv(sendArray.getArray(), recvLens.getArray(),
                recvDispls.getArray(), MPI_INT, receiveArray.getArray(),
                sendLens.getArray(), sendDispls.getArray(), MPI_INT, comm);

  // ###
  // now complete the newMap
  // ###

  for (i = 0; i < totToSend; ++i) {
#ifdef DEBUG_CONTROLLER
    assert(receiveArray[i] >= 0 && receiveArray[i] < numTotalVertices);
    assert(copyOfSendArray[i] >= 0 && copyOfSendArray[i] < numOrigLocVerts);
#endif
    auxArray[copyOfSendArray[i]] = receiveArray[i];
  }

#ifdef DEBUG_CONTROLLER
  for (i = 0; i < numOrigLocVerts; ++i)
    assert(auxArray[i] >= 0 && auxArray[i] < numTotalVertices);
#endif

  mapToOrigVerts.setArray(auxArray, numOrigLocVerts);
}

void ParaVCycleController::resetStructs() {
  hgraph->resetVectors();

  mapToInterVerts.setLength(0);

#ifdef DEBUG_CONTROLLER
  assert(hgraphs.getNumElem() == 0);
  assert(bestVCyclePartition.getNumElem() == 0);
#endif
}

#endif
