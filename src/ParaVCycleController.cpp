
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
#include "data_structures/complete_binary_tree.hpp"

namespace ds = parkway::data_structures;

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

  mapToInterVerts.reserve(0);
  mapToOrigVerts.reserve(0);
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

    array = bestPartVector->data();

    for (i = 0; i < numLocalVertices; ++i)
      array[i] = pVector[i];

    minLocCurrVertId.push(minVertexIndex);
    numLocCurrVerts.push(numLocalVertices);
    bestVCyclePartition.push(bestPartVector);
  } else {
    bestPartVector = bestVCyclePartition.top();

#ifdef DEBUG_CONTROLLER
    assert(bestPartVector);
#endif

    bestPartVector->reserve(numLocalVertices);
    array = bestPartVector->data();

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

  array = bestPartVector->data();

  dynamic_array<int> minStoredIDs(processors_);
  dynamic_array<int> localSendArrayVerts;
  dynamic_array<int> procs(processors_);

  MPI_Allgather(&minStoredVertexIndex, 1, MPI_INT, minStoredIDs.data(), 1,
                MPI_INT, comm);

  for (i = 0; i < processors_; ++i)
    procs[i] = i;

  ds::complete_binary_tree <int> storedVtoProc(procs.data(),
                                               minStoredIDs.data(), processors_);

  for (i = 0; i < processors_; ++i)
    send_lens_[i] = 0;

  for (i = 0; i < numLocalVertices; ++i) {
    origV = myVtoOrigV[i];
#ifdef DEBUG_CONTROLLER
    assert(origV >= 0 && origV < numTotVertices);
#endif

    if (origV >= minStoredVertexIndex && origV < maxStoredVertexIndex) {
      myPartVector[i] = array[origV - minStoredVertexIndex];
    } else {
      j = storedVtoProc.root_value(origV);

      data_out_sets_[j]->assign(send_lens_[j]++, origV);
      data_out_sets_[j]->assign(send_lens_[j]++, i);
    }
  }

  ij = 0;

  for (i = 0; i < processors_; ++i) {
    send_lens_[i] = Shiftr(send_lens_[i], 1);
    send_displs_[i] = ij;
    ij += send_lens_[i];
  }

  send_array_.reserve(ij);
  localSendArrayVerts.reserve(ij);
  totToSend = ij;

  j = 0;

  for (i = 0; i < processors_; ++i) {
    dataOutArray = data_out_sets_[i]->data();
    arrayLen = Shiftl(send_lens_[i], 1);

    for (ij = 0; ij < arrayLen;) {
      send_array_[j] = dataOutArray[ij++];
      localSendArrayVerts[j++] = dataOutArray[ij++];
    }
  }
#ifdef DEBUG_CONTROLLER
  assert(j == totToSend);
#endif

  MPI_Alltoall(send_lens_.data(), 1, MPI_INT, receive_lens_.data(), 1, MPI_INT,
               comm);

  j = 0;
  for (i = 0; i < processors_; ++i) {
    receive_displs_[i] = j;
    j += receive_lens_[i];
  }

  receive_array_.reserve(j);
  totToRecv = j;

  MPI_Alltoallv(send_array_.data(), send_lens_.data(),
                send_displs_.data(), MPI_INT, receive_array_.data(),
                receive_lens_.data(), receive_displs_.data(), MPI_INT, comm);

  // ###
  // now have received all requests and sent out our requests
  // the reply communication will have the dual dimensions
  // ###

  send_array_.reserve(totToRecv);

  for (i = 0; i < totToRecv; ++i) {
#ifdef DEBUG_CONTROLLER
    assert(receive_array_[i] >= minStoredVertexIndex &&
           receive_array_[i] < maxStoredVertexIndex);
#endif

    send_array_[i] = array[receive_array_[i] - minStoredVertexIndex];
  }

  receive_array_.reserve(totToSend);

  MPI_Alltoallv(send_array_.data(), receive_lens_.data(),
                receive_displs_.data(), MPI_INT, receive_array_.data(),
                send_lens_.data(), send_displs_.data(), MPI_INT, comm);

  for (i = 0; i < totToSend; ++i) {
    myPartVector[localSendArrayVerts[i]] = receive_array_[i];
  }

  h.setCut(0, cut);

#ifdef DEBUG_CONTROLLER
  h.checkPartitions(numTotalParts, maxPartWt, comm);
#endif

  DynaMem::deletePtr<IntArray>(bestPartVector);
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

  if (mapToInterVerts.capacity() != 0) {
    fG.setNumberPartitions(1);

    int numLocalCoarseVertices = cG.getNumLocalVertices();
    int numTotalCoarseVertices = cG.getNumTotalVertices();
    int coarseCut = cG.getCut(0);

    fG.setCut(0, coarseCut);

#ifdef DEBUG_CONTROLLER
    assert(numLocalCoarseVertices == mapToInterVerts.capacity());
    if (rank_ != processors_ - 1)
      assert(numLocalCoarseVertices == numTotalCoarseVertices / processors_);
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
    int cVertPerProc = numTotalCoarseVertices / processors_;
    int cVertex;
    int vPart;
    int totToRecv;
    int totToSend;
    int sendLength;

    int j;
    int ij;

    dynamic_array<int> interGraphPVector(numLocalCoarseVertices);
    dynamic_array<int> requestingLocalVerts;

    for (i = 0; i < processors_; ++i)
      send_lens_[i] = 0;

    for (i = 0; i < numLocalCoarseVertices; ++i) {
      cVertex = mapToInterVerts[i];
      vPart = coarsePartVector[i];

#ifdef DEBUG_CONTROLLER
      assert(cVertex >= 0 && cVertex < numTotalCoarseVertices);
      assert(vPart >= 0 && vPart < numTotalParts);
#endif

      ij = min(cVertex / cVertPerProc, processors_ - 1);

      if (ij == rank_) {
        interGraphPVector[cVertex - minCoarseVertexId] = vPart;
      } else {
        data_out_sets_[ij]->assign(send_lens_[ij]++, cVertex);
        data_out_sets_[ij]->assign(send_lens_[ij]++, vPart);
      }
    }

    ij = 0;

    for (i = 0; i < processors_; ++i) {
      send_displs_[i] = ij;
      ij += send_lens_[i];
    }

    send_array_.reserve(ij);
    totToSend = ij;
    ij = 0;

    for (i = 0; i < processors_; ++i) {
      j = 0;
      sendLength = send_lens_[i];
      array = data_out_sets_[i]->data();

      while (j < sendLength) {
        send_array_[ij++] = array[j++];
      }
    }

#ifdef DEBUG_CONTROLLER
    assert(ij == totToSend);
#endif

    // ###
    // get dimension and carry
    // out the communication
    // ###

    MPI_Alltoall(send_lens_.data(), 1, MPI_INT, receive_lens_.data(), 1,
                 MPI_INT, comm);

    ij = 0;
    for (i = 0; i < processors_; ++i) {
      receive_displs_[i] = ij;
      ij += receive_lens_[i];
    }

    receive_array_.reserve(ij);
    totToRecv = ij;

    MPI_Alltoallv(send_array_.data(), send_lens_.data(),
                  send_displs_.data(), MPI_INT, receive_array_.data(),
                  receive_lens_.data(), receive_displs_.data(), MPI_INT, comm);

    // ###
    // now finish the initialisation
    // of the corse hypergraph partitions
    // ###

    for (i = 0; i < totToRecv;) {
      cVertex = receive_array_[i++];
      vPart = receive_array_[i++];

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

    for (i = 0; i < processors_; ++i)
      send_lens_[i] = 0;

    // ###
    // initialise the local partition structures
    // ###

    for (i = 0; i < numLocalFineVertices; ++i) {
#ifdef DEBUG_CONTROLLER
      assert(fineMatVector[i] >= 0 &&
             fineMatVector[i] < numTotalCoarseVertices);
#endif
      cVertex = fineMatVector[i];
      ij = min(cVertex / cVertPerProc, processors_ - 1);

      if (ij == rank_) {
        finePartVector[i] = interGraphPVector[cVertex - minCoarseVertexId];
      } else {
        data_out_sets_[ij]->assign(send_lens_[ij]++, i);
        data_out_sets_[ij]->assign(send_lens_[ij]++, cVertex);
      }
    }

    // ###
    // prepare to send requests for partition vector values
    // ###

    ij = 0;

    for (i = 0; i < processors_; ++i) {
      send_displs_[i] = ij;
      ij += (Shiftr(send_lens_[i], 1));
    }

    send_array_.reserve(ij);
    requestingLocalVerts.reserve(ij);

    numRequestingLocalVerts = ij;
    totToSend = ij;
    ij = 0;

    for (i = 0; i < processors_; ++i) {
      j = 0;
      sendLength = send_lens_[i];
      array = data_out_sets_[i]->data();

      while (j < sendLength) {
        requestingLocalVerts[ij] = array[j++];
        send_array_[ij++] = array[j++];
      }

      send_lens_[i] = Shiftr(sendLength, 1);
    }

#ifdef DEBUG_CONTROLLER
    assert(ij == totToSend);
#endif

    // ###
    // get dimension and carry
    // out the communication
    // ###

    MPI_Alltoall(send_lens_.data(), 1, MPI_INT, receive_lens_.data(), 1,
                 MPI_INT, comm);

    ij = 0;
    for (i = 0; i < processors_; ++i) {
      receive_displs_[i] = ij;
      ij += receive_lens_[i];
    }

    receive_array_.reserve(ij);
    totToRecv = ij;

    MPI_Alltoallv(send_array_.data(), send_lens_.data(),
                  send_displs_.data(), MPI_INT, receive_array_.data(),
                  receive_lens_.data(), receive_displs_.data(), MPI_INT, comm);

    // ###
    // process the requests for
    // local vertex partitions
    // dimensions of cummunication
    // are reversed!
    // ###

    totToSend = totToRecv;
    send_array_.reserve(totToSend);

    for (i = 0; i < totToRecv; ++i) {
#ifdef DEBUG_CONTROLLER
      assert(receive_array_[i] >= minCoarseVertexId &&
             receive_array_[i] < maxCoarseVertexId);
#endif
      cVertex = receive_array_[i] - minCoarseVertexId;
      send_array_[i] = interGraphPVector[cVertex];
    }

    totToRecv = numRequestingLocalVerts;
    receive_array_.reserve(totToRecv);

    // ###
    // now do the dual communication
    // ###

    MPI_Alltoallv(send_array_.data(), receive_lens_.data(),
                  receive_displs_.data(), MPI_INT, receive_array_.data(),
                  send_lens_.data(), send_displs_.data(), MPI_INT, comm);

    // ###
    // finish off initialising the
    // partition vector
    // ###

    ij = 0;

    for (i = 0; i < numRequestingLocalVerts; ++i) {
      cVertex = requestingLocalVerts[i];
#ifdef DEBUG_CONTROLLER
      assert(cVertex >= 0 && cVertex < numLocalFineVertices);
      assert(receive_array_[ij] >= 0 && receive_array_[ij] < numTotalParts);
#endif
      finePartVector[cVertex] = receive_array_[ij++];
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

  mapToInterVerts.reserve(numLocalFineVertices);
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
  int numVBefPerProc = numTotVertices / processors_;

#ifdef DEBUG_CONTROLLER
  assert(minInterVertIndex == minLocVertIdBefShuff);
  assert(numLocVertBefShuff == mapToInterVerts.capacity());
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

  dynamic_array<int> *newToInter = new dynamic_array<int>(numLocVertAftShuff);
  dynamic_array<int> requestingLocalVerts;

  for (i = 0; i < processors_; ++i)
    send_lens_[i] = 0;

  array = newToInter->data();

  for (i = 0; i < numLocVertAftShuff; ++i) {
    j = vToOrigV[i];
    ij = min(j / numVBefPerProc, processors_ - 1);

    if (ij == rank_) {
      array[i] = mapToInterVerts[j - minInterVertIndex];
    } else {
      data_out_sets_[ij]->assign(send_lens_[ij]++, i);
      data_out_sets_[ij]->assign(send_lens_[ij]++, j);
    }
  }

  ij = 0;
  for (i = 0; i < processors_; ++i) {
    send_displs_[i] = ij;
    ij += (Shiftr(send_lens_[i], 1));
  }

  send_array_.reserve(ij);
  requestingLocalVerts.reserve(ij);
  numRequestingLocalVerts = ij;
  totToSend = ij;
  ij = 0;

  for (i = 0; i < processors_; ++i) {
    j = 0;
    sendLength = send_lens_[i];
    array = data_out_sets_[i]->data();

    while (j < sendLength) {
      requestingLocalVerts[ij] = array[j++];
      send_array_[ij++] = array[j++];
    }

    send_lens_[i] = Shiftr(send_lens_[i], 1);
  }

#ifdef DEBUG_CONTROLLER
  assert(ij == totToSend);
#endif

  // ###
  // get dimension and carry
  // out the communication
  // ###

  MPI_Alltoall(send_lens_.data(), 1, MPI_INT, receive_lens_.data(), 1, MPI_INT,
               comm);

  ij = 0;
  for (i = 0; i < processors_; ++i) {
    receive_displs_[i] = ij;
    ij += receive_lens_[i];
  }

  receive_array_.reserve(ij);
  totToRecv = ij;

  MPI_Alltoallv(send_array_.data(), send_lens_.data(),
                send_displs_.data(), MPI_INT, receive_array_.data(),
                receive_lens_.data(), receive_displs_.data(), MPI_INT, comm);

  // ###
  // process the requests for
  // local mapToInterVerts
  // dimensions of cummunication
  // are reversed!
  // ###

  totToSend = totToRecv;
  send_array_.reserve(totToSend);

  for (i = 0; i < totToRecv; ++i) {
#ifdef DEBUG_CONTROLLER
    assert(receive_array_[i] >= minLocVertIdBefShuff &&
           receive_array_[i] < maxLocVertIdBefShuff);
#endif
    j = receive_array_[i] - minLocVertIdBefShuff;
    send_array_[i] = mapToInterVerts[j];
  }

  totToRecv = numRequestingLocalVerts;
  receive_array_.reserve(totToRecv);

  // ###
  // now do the dual communication
  // ###

  MPI_Alltoallv(send_array_.data(), receive_lens_.data(),
                receive_displs_.data(), MPI_INT, receive_array_.data(),
                send_lens_.data(), send_displs_.data(), MPI_INT, comm);

  // ###
  // finish off initialising the
  // new mapToInterVerts vector
  // ###

  array = newToInter->data();

  for (i = 0; i < numRequestingLocalVerts; ++i) {
    j = requestingLocalVerts[i];
#ifdef DEBUG_CONTROLLER
    assert(j >= 0 && j < numLocVertAftShuff);
    assert(receive_array_[i] >= 0 && receive_array_[i] < numTotVertices);
#endif
    array[j] = receive_array_[i];
  }

  minInterVertIndex = minLocVertIdAftShuff;
  mapToInterVerts.set_data(array, numLocVertAftShuff);
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
  //int numVBefPerProc = numTotVertices / processors_;
  int numVPerProc = numTotVertices / processors_;
  int minLocVertexIndex = h.getMinVertexIndex();
  int numLocVertices = h.getNumLocalVertices();

#  ifdef DEBUG_CONTROLLER
  int maxLocVertexIndex = minLocVertexIndex+numLocVertices;
#  endif

#  ifdef DEBUG_CONTROLLER
  assert(numLocVertices == mapToInterVerts.capacity());
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
  int *data_;

  dynamic_array<int>* newToInter = new dynamic_array<int>(numLocVertices);
  dynamic_array<int> requestingLocalVerts;

  for (i=0;i<processors_;++i)
    send_lens_[i] = 0;

  data_ = newToInter->getArray();

  for (i=0;i<numLocVertices;++i)
    {
      j = vToOrigV[i];
      ij = min(j / numVPerProc, processors_-1);

      if(ij == rank_)
        {
          data_[i] = mapToInterVerts[j-minInterVertIndex];
        }
      else
        {
          data_out_sets_[ij]->assign(send_lens_[ij]++, i);
          data_out_sets_[ij]->assign(send_lens_[ij]++, j);
        }
    }

  ij = 0;
  for (i=0;i<processors_;++i)
    {
      send_displs_[i] = ij;
      ij += (Shiftr(send_lens_[i],1));
    }

  send_array_.setLength(ij);
  requestingLocalVerts.setLength(ij);
  numRequestingLocalVerts = ij;
  totToSend = ij;
  ij = 0;

  for (i=0;i<processors_;++i)
    {
      j=0;
      sendLength = send_lens_[i];
      data_ = data_out_sets_[i]->getArray();

      while(j < sendLength)
        {
          requestingLocalVerts[ij] = data_[j++];
          send_array_[ij++] = data_[j++];
        }

      send_lens_[i] = Shiftr(send_lens_[i],1);
    }

#  ifdef DEBUG_CONTROLLER
  assert(ij == totToSend);
#  endif

  // ###
  // get dimension and carry
  // out the communication
  // ###

  MPI_Alltoall(send_lens_.getArray(), 1, MPI_INT, receive_lens_.getArray(), 1, MPI_INT,
comm);

  ij = 0;
  for (i=0;i<processors_;++i)
    {
      receive_displs_[i] = ij;
      ij += receive_lens_[i];
    }

  receive_array_.setLength(ij);
  totToRecv = ij;

  MPI_Alltoallv(send_array_.getArray(), send_lens_.getArray(),
send_displs_.getArray(), MPI_INT, receive_array_.getArray(), receive_lens_.getArray(),
receive_displs_.getArray(), MPI_INT, comm);

  // ###
  // process the requests for
  // local mapToInterVerts
  // dimensions of cummunication
  // are reversed!
  // ###

  totToSend = totToRecv;
  send_array_.setLength(totToSend);

  for (i=0;i<totToRecv;++i)
    {
#  ifdef DEBUG_CONTROLLER
      assert(receive_array_[i] >= minLocVertexIndex && receive_array_[i] <
maxLocVertexIndex);
#  endif
      j = receive_array_[i]-minLocVertexIndex;//minLocVertIdBefShuff;
      send_array_[i] = mapToInterVerts[j];
    }

  totToRecv = numRequestingLocalVerts;
  receive_array_.reserve(totToRecv);

  // ###
  // now do the dual communication
  // ###

  MPI_Alltoallv(send_array_.getArray(), receive_lens_.getArray(),
receive_displs_.getArray(), MPI_INT, receive_array_.data_(), send_lens_.getArray(),
send_displs_.getArray(), MPI_INT, comm);

  // ###
  // finish off initialising the
  // new mapToInterVerts vector
  // ###

  data_ = newToInter->getArray();

  for (i=0;i<numRequestingLocalVerts;++i)
    {
      j = requestingLocalVerts[i];
#  ifdef DEBUG_CONTROLLER
      assert(j >= 0 && j < numLocVertices);
      assert(receive_array_[i] >= 0 && receive_array_[i] < numTotVertices);
#  endif
      data_[j] = receive_array_[i];
    }

  minInterVertIndex = minLocVertexIndex;
  mapToInterVerts.set_data(data_,numLocVertices);
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
  int vPerProc = numTotVertices / processors_;
  int numMyNewVertices;

  if (rank_ != processors_ - 1)
    numMyNewVertices = vPerProc;
  else
    numMyNewVertices = vPerProc + Mod(numTotVertices, processors_);

  dynamic_array<int> minNewIndex(processors_);
  dynamic_array<int> maxNewIndex(processors_);

#ifdef DEBUG_CONTROLLER
  assert(minVertexIndex == minInterVertIndex);
#endif

  h.shiftVerticesToBalance(comm);

  // ###
  // mapToInterVerts = send_array_
  // newToInter = receive_array_
  // ###

  dynamic_array<int> *newToInter = new dynamic_array<int>(numMyNewVertices);

  for (i = 0; i < processors_; ++i) {
    if (i == 0) {
      minNewIndex[i] = 0;
      maxNewIndex[i] = vPerProc;
    } else {
      minNewIndex[i] = maxNewIndex[i - 1];
      if (i == processors_ - 1)
        maxNewIndex[i] = numTotVertices;
      else
        maxNewIndex[i] = minNewIndex[i] + vPerProc;
    }
  }

  for (i = 0; i < processors_; ++i) {
    if (i == 0)
      send_displs_[i] = 0;
    else
      send_displs_[i] = send_displs_[i - 1] + send_lens_[i - 1];

    send_lens_[i] =
        max(numLocalVertices - (max(maxVertexIndex - maxNewIndex[i], 0) +
                                max(minNewIndex[i] - minVertexIndex, 0)),
            0);
  }

  MPI_Alltoall(send_lens_.data(), 1, MPI_INT, receive_lens_.data(), 1, MPI_INT,
               comm);

  j = 0;
  for (i = 0; i < processors_; ++i) {
    receive_displs_[i] = j;
    j += receive_lens_[i];
  }
#ifdef DEBUG_CONTROLLER
  assert(j == numMyNewVertices);
#endif

  MPI_Alltoallv(mapToInterVerts.data(), send_lens_.data(),
                send_displs_.data(), MPI_INT, newToInter->data(),
                receive_lens_.data(), receive_displs_.data(), MPI_INT, comm);

  minInterVertIndex = minNewIndex[rank_];
  mapToInterVerts.set_data(newToInter->data(), numMyNewVertices);
}

void ParaVCycleController::updateMapToOrigVerts(MPI_Comm comm) {
#ifdef DEBUG_CONTROLLER
  int numLocalVertices = hgraph->getNumLocalVertices();
#endif

  int minLocVertIndex = hgraph->getMinVertexIndex();
  int numTotalVertices = hgraph->getNumTotalVertices();
  int vertPerProc = numTotalVertices / processors_;
  int totToRecv;
  int totToSend;
  int sendLength;
  int vertex;

  int i;
  int j;
  int ij;

#ifdef DEBUG_CONTROLLER
  assert(mapToOrigVerts.getLength() == numOrigLocVerts);
  assert(mapToInterVerts.capacity() == numOrigLocVerts);
  if (rank_ != processors_ - 1)
    assert(numLocalVertices == vertPerProc);
#endif

  dynamic_array<int> *newMapToOrig = new dynamic_array<int>(numOrigLocVerts);
  dynamic_array<int> copyOfSendArray;

  int *auxArray = newMapToOrig->data();
  int *array;

  for (i = 0; i < processors_; ++i)
    send_lens_[i] = 0;

  for (i = 0; i < numOrigLocVerts; i++) {
    vertex = mapToInterVerts[i];

#ifdef DEBUG_CONTROLLER
    assert(vertex >= 0 && vertex < numTotalVertices);
#endif

    ij = min(vertex / vertPerProc, processors_ - 1);

    if (ij == rank_) {
      auxArray[i] = mapToOrigVerts[vertex - minLocVertIndex];
    } else {
      data_out_sets_[ij]->assign(send_lens_[ij]++, i);
      data_out_sets_[ij]->assign(send_lens_[ij]++, vertex);
    }
  }

  ij = 0;

  for (i = 0; i < processors_; ++i) {
    send_displs_[i] = ij;
    ij += Shiftr(send_lens_[i], 1);
  }

  send_array_.reserve(ij);
  copyOfSendArray.reserve(ij);
  totToSend = ij;
  ij = 0;

  for (i = 0; i < processors_; ++i) {
    j = 0;
    sendLength = send_lens_[i];
    array = data_out_sets_[i]->data();

    while (j < sendLength) {
      copyOfSendArray[ij] = array[j++];
      send_array_[ij++] = array[j++];
    }

    send_lens_[i] = Shiftr(send_lens_[i], 1);
  }

#ifdef DEBUG_CONTROLLER
  assert(ij == totToSend);
#endif

  // ###
  // get dimension and carry
  // out the communication
  // ###

  MPI_Alltoall(send_lens_.data(), 1, MPI_INT, receive_lens_.data(), 1, MPI_INT,
               comm);

  ij = 0;
  for (i = 0; i < processors_; ++i) {
    receive_displs_[i] = ij;
    ij += receive_lens_[i];
  }

  receive_array_.reserve(ij);
  totToRecv = ij;

  MPI_Alltoallv(send_array_.data(), send_lens_.data(),
                send_displs_.data(), MPI_INT, receive_array_.data(),
                receive_lens_.data(), receive_displs_.data(), MPI_INT, comm);

  // ###
  // now initialise the bestPartition data_
  // using the data_ in receive_array_
  // ###

  send_array_.reserve(totToRecv);

  for (i = 0; i < totToRecv; ++i) {
#ifdef DEBUG_CONTROLLER
    assert(receive_array_[i] >= minLocVertIndex &&
           receive_array_[i] < minLocVertIndex + numOrigLocVerts);
#endif
    send_array_[i] = mapToOrigVerts[receive_array_[i] - minLocVertIndex];
  }

  receive_array_.reserve(totToSend);

  // ###
  // now do 'dual' communication
  // ###

  MPI_Alltoallv(send_array_.data(), receive_lens_.data(),
                receive_displs_.data(), MPI_INT, receive_array_.data(),
                send_lens_.data(), send_displs_.data(), MPI_INT, comm);

  // ###
  // now complete the newMap
  // ###

  for (i = 0; i < totToSend; ++i) {
#ifdef DEBUG_CONTROLLER
    assert(receive_array_[i] >= 0 && receive_array_[i] < numTotalVertices);
    assert(copyOfSendArray[i] >= 0 && copyOfSendArray[i] < numOrigLocVerts);
#endif
    auxArray[copyOfSendArray[i]] = receive_array_[i];
  }

#ifdef DEBUG_CONTROLLER
  for (i = 0; i < numOrigLocVerts; ++i)
    assert(auxArray[i] >= 0 && auxArray[i] < numTotalVertices);
#endif

  mapToOrigVerts.set_data(auxArray, numOrigLocVerts);
}

void ParaVCycleController::resetStructs() {
  hgraph->resetVectors();

  mapToInterVerts.reserve(0);

#ifdef DEBUG_CONTROLLER
  assert(hgraphs.getNumElem() == 0);
  assert(bestVCyclePartition.getNumElem() == 0);
#endif
}

#endif
