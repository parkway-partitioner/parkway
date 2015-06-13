
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

ParaController::ParaController(parallel_coarsener &c, ParaRefiner &r,
                               SeqController &con, int rank, int nP,
                               int percentile, int inc, int approxRef,
                               ostream &out)
    : global_communicator(rank, nP), out_stream(out), coarsener(c), refiner(r),
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

  bestPartition.reserve(0);
}

ParaController::~ParaController() {}

/*
void ParaController::setShuffleFile(const char *filename)
{
  shuffleFile.reserve(strlen(filename)+1);
  strcpy(shuffleFile.data_(),filename);
}
*/

void ParaController::initMapToOrigVerts() {
  int i;
  int j = hgraph->minimum_vertex_index();

  mapToOrigVerts.reserve(numOrigLocVerts);

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
    int numVPerProc = hgraph->total_number_of_vertices() / processors_;
    int myOffset = rank_ * numVPerProc;

    char message[512];
    ifstream in_stream;

    in_stream.open(filename, ifstream::in | ifstream::binary);

    if (!in_stream.is_open()) {
      sprintf(message, "p[%d] could not open partition file %s\n", rank_,
              filename);
      out_stream << message;
      MPI_Abort(comm, 0);
    }

    shufflePartition.reserve(numOrigLocVerts);

    in_stream.seekg(myOffset * sizeof(int), ifstream::beg);
    len = numOrigLocVerts * sizeof(int);
    in_stream.read((char *)(shufflePartition.data()), len);

    if (in_stream.gcount() != len) {
      sprintf(message, "p[%d] could not read in %d elements\n", rank_,
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

  int minLocVertIndex = hgraph->minimum_vertex_index();
  int numTotalVertices = hgraph->total_number_of_vertices();
  int vertPerProc = numTotalVertices / processors_;
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
  if (rank_ != processors_ - 1)
    assert(numLocalVertices == vertPerProc);
#endif

  int *mapToHgraphVerts = mapToOrigVerts.data();
  int *auxArray;

  for (i = 0; i < processors_; ++i)
    send_lens_[i] = 0;

  bestPartition.reserve(numV);

  for (i = 0; i < numV; ++i) {
    vertex = mapToHgraphVerts[i];
    vPart = array[i];

#ifdef DEBUG_CONTROLLER
    assert(vertex >= 0 && vertex < numTotalVertices);
    assert(vPart >= 0 && vPart < numTotalParts);
#endif

    ij = min(vertex / vertPerProc, processors_ - 1);
    assert(ij < processors_);
    if (ij == rank_) {
      bestPartition[vertex - minLocVertIndex] = vPart;
    } else {
      data_out_sets_[ij]->assign(send_lens_[ij]++, vertex);
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
    auxArray = data_out_sets_[i]->data();

    while (j < sendLength) {
      send_array_[ij++] = auxArray[j++];
    }
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

  for (i = 0; i < totToRecv;) {
    vertex = receive_array_[i++];
    vPart = receive_array_[i++];

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

  if (rank_ == 0)
    remove(filename);

  MPI_Barrier(comm);

  for (i = 0; i < processors_; ++i) {
    if (rank_ == i) {
      out.open(filename, ofstream::out | ofstream::app | ofstream::binary);

      if (!out.is_open()) {
        sprintf(message, "p[%d] cannot open %s\n", rank_, filename);
        out_stream << message;
      } else
        out.write((char *)(bestPartition.data()),
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

  locGraphWt = hgraph->vertex_weight();

  MPI_Allreduce(&locGraphWt, &totGraphWt, 1, MPI_INT, MPI_SUM, comm);

  avePartWt = static_cast<double>(totGraphWt) / numTotalParts;
  maxPartWt = static_cast<int>(floor(avePartWt + avePartWt * balConstraint));
  maxVertWt = static_cast<int>(floor(avePartWt * balConstraint));

  coarsener.set_maximum_vertex_weight(maxVertWt);
  coarsener.set_total_hypergraph_weight(totGraphWt);
  seqController.setMaxVertexWt(maxVertWt);
}

#ifdef DEBUG_TABLES
void ParaController::printHashMemUse() {
  // write_log(rank_, "HedgeIndex diff: %d",
  // HedgeIndexEntry::getNumAllocated()-HedgeIndexEntry::getNumDeleted());
  // write_log(rank_, "MatchEntry diff: %d",
  // MatchEntry::getNumAllocated()-MatchEntry::getNumDeleted());
  write_log(rank_, "MatchRequestEntry diff: %d",
            MatchRequestEntry::getNumAllocated() -
                MatchRequestEntry::getNumDeleted());
  // write_log(rank_, "IndexEntry diff: %d",
  // IndexEntry::getNumAllocated()-IndexEntry::getNumDeleted());
  // write_log(rank_, "ConnVertData diff: %d",
  // ConnVertData::getNumAllocated()-ConnVertData::getNumDeleted());
  // write_log(rank_, "VertexPartEntry diff: %d",
  // VertexPartEntry::getNumAllocated()-VertexPartEntry::getNumDeleted());
  // write_log(rank_, "DuplRemEntry diff: %d",
  // DuplRemEntry::getNumAllocated()-DuplRemEntry::getNumDeleted());
  // write_log(rank_, "PVectorEntry diff: %d",
  // PVectorEntry::getNumAllocated()-PVectorEntry::getNumDeleted());
}
#endif

#endif
