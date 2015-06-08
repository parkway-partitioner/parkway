
#ifndef _PARA_CONTROLLER_CPP
#define _PARA_CONTROLLER_CPP

// ### ParaController.cpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 4/1/2005: Last Modified
//
// ###

#include "ParaController.hpp"

ParaController::ParaController(ParaCoarsener &c, ParaRefiner &r,
                               SeqController &con, int rank, int nP,
                               int percentile, int inc, int approxRef,
                               ostream &out)
    : GlobalCommunicator(rank, nP), out_stream(out), coarsener(c), refiner(r),
      seqController(con) {
  shuffled = 0;
  numParaRuns = -1;
  numTotalParts = -1;
  dispOption = -1;
  writePartitionToFile = -1;

  startPercentile = percentile;
  percentileIncrement = inc;
  approxRefine = approxRef;

  maxPartWt = 0;
  numOrigLocVerts = 0;
  keepPartitionsWithin = 0;
  reductionInKeepThreshold = 1.0;
  balConstraint = 0;
  bestCutsize = 0;
  worstCutsize = 0;
  aveCutSize = 0;
  startTime = 0;
  totCoaTime = 0;
  totSeqTime = 0;
  totRefTime = 0;
  totalTime = 0;
  accumulator = 0;

  hgraph = nullptr;

  bestPartition.setLength(0);
}

ParaController::~ParaController() {}

/*
void ParaController::setShuffleFile(const char *filename)
{
  shuffleFile.setLength(strlen(filename)+1);
  strcpy(shuffleFile.getArray(),filename);
}
*/

void ParaController::initMapToOrigVerts() {
  int i;
  int j = hgraph->getMinVertexIndex();

  mapToOrigVerts.setLength(numOrigLocVerts);

  for (i = 0; i < numOrigLocVerts; ++i)
    mapToOrigVerts[i] = j + i;

  /*
  if(shuffled)
    {
      int *mapToOrig = hgraph->getToOrigVArray();

      for (i=0;i<numOrigLocVerts;++i)
        mapToOrigVerts[i] = mapToOrig[i];
    }
  else
    {
  */
}

void ParaController::setPrescribedPartition(const char *filename,
                                            MPI_Comm comm) {
#ifdef DEBUG_CONTROLLER
  assert(hgraph);
#endif

  if (shuffled == 2) {
    int len;
    int numVPerProc = hgraph->getNumTotalVertices() / numProcs;
    int myOffset = myRank * numVPerProc;

    char message[512];
    ifstream in_stream;

    in_stream.open(filename, ifstream::in | ifstream::binary);

    if (!in_stream.is_open()) {
      sprintf(message, "p[%d] could not open partition file %s\n", myRank,
              filename);
      out_stream << message;
      MPI_Abort(comm, 0);
    }

    shufflePartition.setLength(numOrigLocVerts);

    in_stream.seekg(myOffset * sizeof(int), ifstream::beg);
    len = numOrigLocVerts * sizeof(int);
    in_stream.read((char *)(shufflePartition.getArray()), len);

    if (in_stream.gcount() != len) {
      sprintf(message, "p[%d] could not read in %d elements\n", myRank,
              numOrigLocVerts);
      out_stream << message;
      MPI_Abort(comm, 0);
    }

    in_stream.close();
  }
}

void ParaController::storeBestPartition(int numV, const int *array,
                                        MPI_Comm comm) {

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
  int vPart;

  int i;
  int j;
  int ij;

#ifdef DEBUG_CONTROLLER
  assert(numLocalVertices == numV);
  assert(mapToOrigVerts.getLength() == numV);
  if (myRank != numProcs - 1)
    assert(numLocalVertices == vertPerProc);
#endif

  int *mapToHgraphVerts = mapToOrigVerts.getArray();
  int *auxArray;

  for (i = 0; i < numProcs; ++i)
    sendLens[i] = 0;

  bestPartition.setLength(numV);

  for (i = 0; i < numV; ++i) {
    vertex = mapToHgraphVerts[i];
    vPart = array[i];

#ifdef DEBUG_CONTROLLER
    assert(vertex >= 0 && vertex < numTotalVertices);
    assert(vPart >= 0 && vPart < numTotalParts);
#endif

    ij = min(vertex / vertPerProc, numProcs - 1);
    assert(ij < numProcs);
    if (ij == myRank) {
      bestPartition[vertex - minLocVertIndex] = vPart;
    } else {
      dataOutSets[ij]->assign(sendLens[ij]++, vertex);
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
    auxArray = dataOutSets[i]->getArray();

    while (j < sendLength) {
      sendArray[ij++] = auxArray[j++];
    }
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

  for (i = 0; i < totToRecv;) {
    vertex = receiveArray[i++];
    vPart = receiveArray[i++];

#ifdef DEBUG_CONTROLLER
    assert(vertex >= minLocVertIndex && vertex < minLocVertIndex + numV);
    assert(vPart >= 0 && vPart < numTotalParts);
#endif

    bestPartition[vertex - minLocVertIndex] = vPart;
  }

#ifdef DEBUG_CONTROLLER
  for (i = 0; i < numV; ++i)
    assert(bestPartition[i] >= 0 && bestPartition[i] < numTotalParts);
#endif
}

void ParaController::partitionToFile(const char *filename,
                                     MPI_Comm comm) const {
  int i;

  char message[512];
  ofstream out;

  if (myRank == 0)
    remove(filename);

  MPI_Barrier(comm);

  for (i = 0; i < numProcs; ++i) {
    if (myRank == i) {
      out.open(filename, ofstream::out | ofstream::app | ofstream::binary);

      if (!out.is_open()) {
        sprintf(message, "p[%d] cannot open %s\n", myRank, filename);
        out_stream << message;
      } else
        out.write((char *)(bestPartition.getArray()),
                  sizeof(int) * numOrigLocVerts);

      out.close();
    }
    MPI_Barrier(comm);
  }
}

void ParaController::copyOutPartition(int numVertices,
                                      int *pVector) const {
  int i;

  for (i = 0; i < numVertices; ++i)
    pVector[i] = bestPartition[i];
}

void ParaController::setWeightConstraints(MPI_Comm comm) {
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
  coarsener.setTotGraphWt(totGraphWt);
  seqController.setMaxVertexWt(maxVertWt);
}

#ifdef DEBUG_TABLES
void ParaController::printHashMemUse() {
  // write_log(myRank, "HedgeIndex diff: %d",
  // HedgeIndexEntry::getNumAllocated()-HedgeIndexEntry::getNumDeleted());
  // write_log(myRank, "MatchEntry diff: %d",
  // MatchEntry::getNumAllocated()-MatchEntry::getNumDeleted());
  write_log(myRank, "MatchRequestEntry diff: %d",
            MatchRequestEntry::getNumAllocated() -
                MatchRequestEntry::getNumDeleted());
  // write_log(myRank, "IndexEntry diff: %d",
  // IndexEntry::getNumAllocated()-IndexEntry::getNumDeleted());
  // write_log(myRank, "ConnVertData diff: %d",
  // ConnVertData::getNumAllocated()-ConnVertData::getNumDeleted());
  // write_log(myRank, "VertexPartEntry diff: %d",
  // VertexPartEntry::getNumAllocated()-VertexPartEntry::getNumDeleted());
  // write_log(myRank, "DuplRemEntry diff: %d",
  // DuplRemEntry::getNumAllocated()-DuplRemEntry::getNumDeleted());
  // write_log(myRank, "PVectorEntry diff: %d",
  // PVectorEntry::getNumAllocated()-PVectorEntry::getNumDeleted());
}
#endif

#endif
