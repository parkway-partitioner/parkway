#ifndef _PARA_HYPERGRAPH_CPP
#define _PARA_HYPERGRAPH_CPP

// ### ParaHypergraph.cpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 13/4/2004: modified duplicate hyperedge removal
//            algorithm. Now uses hash key to compute
//            destination processor for hyperedge
//            where decision is made on which processor
//            to keep hyperedge.
//
// 17/4/2004: complete reconstruction of data_ structures
//            removed HedgeTable and replaced with a
//            local pinlist and arrays for other
//            hyperedge attributes. Access to hyperedges
//            via hashkey remains and done via hEdgeIndexTable
//
// 3/2/2005: Last Modified
//
// ###

#include "ParaHypergraph.hpp"
#include "data_structures/bit_field.hpp"
#include "data_structures/complete_binary_tree.hpp"
#include "data_structures/map_from_pos_int.hpp"
#include "data_structures/new_hyperedge_index_table.hpp"
#include "Log.h"

namespace ds = parkway::data_structures;

ParaHypergraph::ParaHypergraph(int rank_, int nProcs, int _numLocVerts,
                               int _totVerts, int _minVertIndex, int coarsen,
                               int *wtArray)
    : global_communicator(rank_, nProcs),
      doNotCoarsen(coarsen),
      numTotalVertices(_totVerts),
      numLocalVertices(_numLocVerts),
      minVertexIndex(_minVertIndex),
      localVertexWt(0),
      numPartitions(0) {

  vWeight.set_data(wtArray, numLocalVertices);
  matchVector.reserve(numLocalVertices);
  vToOrigV.reserve(0);

  for (int i = 0; i < numLocalVertices; ++i) {
    matchVector[i] = -1;
    localVertexWt += vWeight[i];
  }
}

ParaHypergraph::ParaHypergraph(int rank_, int nProcs, int _numLocVerts,
                               int _totVerts, int _minVertIndex, int coarsen,
                               int cut, int *wtArray, int *partArray)
    : global_communicator(rank_, nProcs),
      doNotCoarsen(coarsen),
      numTotalVertices(_totVerts),
      numLocalVertices(_numLocVerts),
      minVertexIndex(_minVertIndex),
      localVertexWt(0),
      numPartitions(1) {

  vToOrigV.reserve(0);
  vWeight.set_data(wtArray, numLocalVertices);
  partitionVector.set_data(partArray, numLocalVertices);
  matchVector.reserve(numLocalVertices);
  partitionOffsetsVector.reserve(numPartitions + 1);
  partitionCutsizesVector.reserve(numPartitions);

  partitionCutsizesVector[0] = cut;
  partitionOffsetsVector[0] = 0;
  partitionOffsetsVector[1] = numLocalVertices;

  for (int i = 0; i < numLocalVertices; ++i) {
    matchVector[i] = -1;
    localVertexWt += vWeight[i];
  }
}

ParaHypergraph::ParaHypergraph(int rank_, int nProcs, const char *filename,
                               int dispOption, std::ostream &out, MPI_Comm comm)
    : global_communicator(rank_, nProcs) {
  hypergraphFromFile(filename, dispOption, out, comm);
}

ParaHypergraph::ParaHypergraph(int rank_, int nProcs, int numLocVerts,
                               int numLocHedges, int maxHedgeLen,
                               const int *vWeights, const int *hEdgeWts,
                               const int *locPinList, const int *offsets,
                               int dispOption, std::ostream &out, MPI_Comm comm)
    : global_communicator(rank_, nProcs) {

  MPI_Allreduce(&numLocVerts, &numTotalVertices, 1, MPI_INT, MPI_SUM, comm);
  MPI_Scan(&numLocVerts, &minVertexIndex, 1, MPI_INT, MPI_SUM, comm);

  numLocalVertices = numLocVerts;
  minVertexIndex -= numLocalVertices;
  localVertexWt = 0;

  vWeight.reserve(numLocalVertices);
  vToOrigV.reserve(0);
  matchVector.reserve(numLocalVertices);

  for (int i = 0; i < numLocalVertices; ++i) {
    vWeight[i] = vWeights[i];
    localVertexWt += vWeight[i];
    matchVector[i] = -1;
  }

  /* hyperedges stored in pinlist */

  int hEdgeIndex = 0;
  int pinCounter = 0;

  int startOffset;
  int endOffset;
  int hEdgeLength;
  for (int i = 0; i < numLocHedges; ++i) {
    startOffset = offsets[i];
    endOffset = offsets[i + 1];
    hEdgeLength = endOffset - startOffset;

    if (hEdgeLength > 1) {
      hEdgeWeights.assign(hEdgeIndex, hEdgeWts[i]);
      hEdgeOffsets.assign(hEdgeIndex++, pinCounter);

      for (int j = startOffset; j < endOffset; ++j) {
#ifdef DEBUG_HYPERGRAPH
        assert(locPinList[j] >= 0 && locPinList[j] < numTotalVertices);
#endif
        localPins.assign(pinCounter++, locPinList[j]);
      }
    }
  }

  hEdgeOffsets.assign(hEdgeIndex, pinCounter);

  numLocalPins = pinCounter;
  numLocalHedges = hEdgeIndex;
  doNotCoarsen = 0;
  numPartitions = 0;

  int i;
  int j;

  if (dispOption > 0) {
    MPI_Reduce(&numLocalPins, &i, 1, MPI_INT, MPI_SUM, 0, comm);
    MPI_Reduce(&numLocalHedges, &j, 1, MPI_INT, MPI_SUM, 0, comm);

    if (rank_ == 0) {
      out << "|--- Hypergraph (as loaded):" << std::endl;
      out << "| |V| = " << numTotalVertices;
      out << " |E| = " << j;
      out << " |Pins| = " << i << std::endl;
      out << "|" << std::endl;
    }
  }

#ifdef MEM_OPT
  hEdgeOffsets.reserve(numLocalHedges + 1);
  hEdgeWeights.reserve(numLocalHedges);
  localPins.reserve(numLocalPins);
#endif
}

ParaHypergraph::~ParaHypergraph() {}

void ParaHypergraph::hypergraphFromFile(const char *filename, int dispOption,
                                        std::ostream &out, MPI_Comm comm) {
  int buffer[3];
  int numHedgesInFile;
  int hEdgeDataLength;
  int hEdgeLength;
  int hEdgeChunk;
  int index;
  int maxHedgeLen;

  int i;
  int j;
  int hEdgeIndex;
  int pinCounter;

  char my_file[512];
  char message[512];

  std::ifstream in_stream;

  dynamic_array<int> hEdgeData;

  sprintf(my_file, "%s-%d", filename, rank_);

  in_stream.open(my_file, std::ifstream::in | std::ifstream::binary);

  if (!in_stream.is_open()) {
    sprintf(message, "p[%d] could not open file %s\n", rank_, my_file);
    out << message;
    MPI_Abort(comm, 0);
  }

  i = sizeof(int) * 3;
  in_stream.read((char *)(&buffer[0]), i);
  if (in_stream.gcount() != i) {
    sprintf(message, "p[%d] could not read in int buffer[3]\n", rank_);
    out << message;
    in_stream.close();
    MPI_Abort(comm, 0);
  }

  numTotalVertices = buffer[0];
  numLocalVertices = buffer[1];
  hEdgeDataLength = buffer[2];

  minVertexIndex = (numTotalVertices / processors_) * rank_;

  vToOrigV.reserve(0);
  vWeight.reserve(numLocalVertices);
  matchVector.reserve(numLocalVertices);
  hEdgeData.reserve(hEdgeDataLength);

  i = sizeof(int) * numLocalVertices;
  in_stream.read((char *)(vWeight.data()), i);
  if (in_stream.gcount() != i) {
    sprintf(message, "p[%d] could not read in %d vertex elements\n", rank_,
            numLocalVertices);
    out << message;
    in_stream.close();
    MPI_Abort(comm, 0);
  }

  // write_log(rank_, "read in vertex weights");

  i = sizeof(int) * hEdgeDataLength;
  in_stream.read((char *)(hEdgeData.data()), i);
  if (in_stream.gcount() != i) {
    sprintf(message, "p[%d] could not read in %d hyperedge data_ elements\n",
            rank_, hEdgeDataLength);
    out << message;
    in_stream.close();
    MPI_Abort(comm, 0);
  }

  in_stream.close();

  numHedgesInFile = 0;
  j = 0;

  for (i = 0; i < hEdgeDataLength; i += hEdgeData[i]) {
    hEdgeIndex = hEdgeData[i] - 2;

    if (hEdgeIndex > j)
      j = hEdgeIndex;

    ++numHedgesInFile;
  }

  MPI_Allreduce(&j, &maxHedgeLen, 1, MPI_INT, MPI_MAX, comm);
  MPI_Reduce(&numHedgesInFile, &i, 1, MPI_INT, MPI_SUM, 0, comm);

  Funct::setMaxHedgeLen(maxHedgeLen);

  if (rank_ == 0 && dispOption > 0) {
    out << "|--- Hypergraph " << filename << " (on file):" << std::endl
        << "| |V| = " << numTotalVertices << " |E| = " << i << std::endl;
  }

  localVertexWt = 0;
  for (i = 0; i < numLocalVertices; ++i) {
    localVertexWt += vWeight[i];
    matchVector[i] = -1;
  }

  // ### hyperedges stored on file in blocks:
  // [0]  =  block capacity_
  // [1]  =  hyperedge weight
  // [2]-[block capacity_-1] = hyperedge vertices
  // ###

  i = 0;
  index = 0;
  hEdgeIndex = 0;
  pinCounter = 0;

  while (i < hEdgeDataLength) {
    hEdgeChunk = hEdgeData[i];
    hEdgeLength = hEdgeChunk - 2;

    if (hEdgeLength > 1) {
      hEdgeWeights.assign(hEdgeIndex, hEdgeData[i + 1]);
      hEdgeOffsets.assign(hEdgeIndex++, pinCounter);

      for (j = 2; j < hEdgeChunk; ++j) {
#ifdef DEBUG_HYPERGRAPH
        assert(hEdgeData[i + j] >= 0 && hEdgeData[i + j] < numTotalVertices);
#endif
        localPins.assign(pinCounter++, hEdgeData[i + j]);
      }
    }

    i += hEdgeChunk;
  }

  hEdgeOffsets.assign(hEdgeIndex, pinCounter);

  numLocalPins = pinCounter;
  numLocalHedges = hEdgeIndex;
  doNotCoarsen = 0;
  numPartitions = 0;

#ifdef DEBUG_HYPERGRAPH
  for (i = 0; i < numLocalHedges; ++i)
    assert(hEdgeWeights[i] > 0);
#endif

  if (dispOption > 0) {
    MPI_Reduce(&numLocalPins, &i, 1, MPI_INT, MPI_SUM, 0, comm);
    MPI_Reduce(&numLocalHedges, &j, 1, MPI_INT, MPI_SUM, 0, comm);

    if (rank_ == 0) {
      out << "|--- Hypergraph " << filename << " (as loaded):" << std::endl
          << "| |V| = " << numTotalVertices << " |E| = " << j
          << " |Pins| = " << i << std::endl
          << "| # Processors = " << processors_ << std::endl
          << "| " << std::endl;
    }
  }

#ifdef MEM_OPT
  hEdgeOffsets.reserve(numLocalHedges + 1);
  hEdgeWeights.reserve(numLocalHedges);
  localPins.reserve(numLocalPins);
#endif
}

void ParaHypergraph::initPartitionFromFile(const char *filename, int numParts,
                                           std::ostream &out, MPI_Comm comm) {
  int i;
  int myOffset;
  int numVPerProc = numTotalVertices / processors_;

  char message[512];

  std::ifstream in_stream;

  numPartitions = 1;

  partitionVector.reserve(numLocalVertices);
  partitionOffsetsVector.reserve(2);
  partitionCutsizesVector.reserve(1);

  partitionOffsetsVector[0] = 0;
  partitionOffsetsVector[1] = numLocalVertices;

  myOffset = 0;
  for (i = 0; i < rank_; ++i)
    myOffset += numVPerProc;

  in_stream.open(filename, std::ifstream::in | std::ifstream::binary);

  if (!in_stream.is_open()) {
    sprintf(message, "p[%d] could not open partition file %s\n", rank_,
            filename);
    out << message;
    MPI_Abort(comm, 0);
  }

  in_stream.seekg(myOffset * sizeof(int), std::ifstream::beg);

  i = numLocalVertices * sizeof(int);

  in_stream.read((char *)(partitionVector.data()), i);
  if (in_stream.gcount() != i) {
    sprintf(message, "p[%d] could not read in %d elements\n", rank_,
            numLocalVertices);
    out << message;
    MPI_Abort(comm, 0);
  }

  in_stream.close();

  partitionCutsizesVector[0] = calcCutsize(numParts, 0, comm);
}

void ParaHypergraph::allocHedgeMem(int numHedges, int numLocPins) {
#ifdef DEBUG_HYPERGRAPH
  assert(numLocalHedges == numHedges);
  assert(numLocalPins == numLocPins);
#endif

  hEdgeOffsets.reserve(numHedges + 1);
  hEdgeWeights.reserve(numHedges);
  localPins.reserve(numLocPins);
}

void ParaHypergraph::contractHyperedges(ParaHypergraph &coarse, MPI_Comm comm) {
  int vFinePerProc;
  int maxLocalVertex = minVertexIndex + numLocalVertices;
  int totalToSend;
  int totalToRecv;
  int arrayLen;
  int p;
  int vertex;
  int startOffset;
  int endOffset;
  int contractedPinListLen;
  int contractedHedgeLen;
  int numContractedHedges;
  int maxLocalCoarseHedgeLen = 0;
  int maxCoarseHedgeLen;

  int *array;

  int i;
  int j;

#ifdef DEBUG_HYPERGRAPH
  int numTotCoarseVerts = coarse.getNumTotalVertices();
  for (i = 0; i < numLocalVertices; ++i)
    assert(matchVector[i] >= 0 && matchVector[i] < numTotCoarseVerts);
#endif

  vFinePerProc = numTotalVertices / processors_;

  dynamic_array<int> contractedPinList;
  dynamic_array<int> contractedHedgeOffsets;
  dynamic_array<int> contractedHedgeWts;
  dynamic_array<int> origContractedPinList(numLocalPins);
  dynamic_array<int> copyOfReq;

  for (i = 0; i < processors_; ++i)
    send_lens_[i] = 0;

  ds::bit_field sentRequests(numTotalVertices);
  sentRequests.unset();

  // ###
  // compute all the requests for remote vertex matches
  // ###

  for (i = 0; i < numLocalPins; ++i) {
    vertex = localPins[i];

    if (vertex < minVertexIndex || vertex >= maxLocalVertex) {
      if (!sentRequests(vertex)) {
        p = std::min(vertex / vFinePerProc, processors_ - 1);
        data_out_sets_[p]->assign(send_lens_[p]++, vertex);
        sentRequests.set(vertex);
      }

      origContractedPinList[i] = -1;
    } else {
      origContractedPinList[i] = matchVector[vertex - minVertexIndex];
    }
  }

  // ###
  // compute number of elements to send to other procs
  // ###

  j = 0;
  for (i = 0; i < processors_; ++i) {
    send_displs_[i] = j;
    j += send_lens_[i];
  }

  send_array_.reserve(j);
  copyOfReq.reserve(j);
  totalToSend = j;

  j = 0;
  for (p = 0; p < processors_; ++p) {
    array = data_out_sets_[p]->data();
    arrayLen = send_lens_[p];

    for (i = 0; i < arrayLen; ++i) {
      vertex = array[i];
      send_array_[j] = vertex;
      copyOfReq[j++] = vertex;
    }
  }

  MPI_Alltoall(send_lens_.data(), 1, MPI_INT, receive_lens_.data(), 1, MPI_INT,
               comm);

  // ###
  // compute number of elements to receive from other procs
  // ###

  j = 0;
  for (i = 0; i < processors_; ++i) {
    receive_displs_[i] = j;
    j += receive_lens_[i];
  }

  receive_array_.reserve(j);
  totalToRecv = j;

  MPI_Alltoallv(send_array_.data(), send_lens_.data(),
                send_displs_.data(), MPI_INT, receive_array_.data(),
                receive_lens_.data(), receive_displs_.data(), MPI_INT, comm);

  // ###
  // now have received all requests and sent out our requests
  // the reply communication will have the dual dimensions
  // ###

  send_array_.reserve(totalToRecv);

  for (i = 0; i < totalToRecv; ++i) {
#ifdef DEBUG_HYPERGRAPH
    assert(receive_array_[i] >= minVertexIndex &&
           receive_array_[i] < maxLocalVertex);
#endif
    send_array_[i] = matchVector[receive_array_[i] - minVertexIndex];
  }

  receive_array_.reserve(totalToSend);

  MPI_Alltoallv(send_array_.data(), receive_lens_.data(),
                receive_displs_.data(), MPI_INT, receive_array_.data(),
                send_lens_.data(), send_displs_.data(), MPI_INT, comm);

  // ###
  // now the requested vertices are in the copyOfReq data_
  // while their corresponding matchVector values are in
  // the corresponding location in the receive_array_
  // ###

  /* depending on number of local pins, choose format for
     storing non-local vertices */

  if (numLocalPins < numTotalVertices / 2) {
    ds::map_from_pos_int<int> storedRequests(numLocalPins);

    for (i = 0; i < totalToSend; ++i)
      storedRequests.insert(copyOfReq[i], receive_array_[i]);

    /* contract remaining local pins */

    for (i = 0; i < numLocalPins; ++i)
      if (origContractedPinList[i] == -1)
        origContractedPinList[i] = storedRequests.get(localPins[i]);
  } else {
    dynamic_array<int> nonLocalMatches(numTotalVertices);

    for (i = 0; i < totalToSend; ++i)
      nonLocalMatches[copyOfReq[i]] = receive_array_[i];

    /* contract remaining local pins */

    for (i = 0; i < numLocalPins; ++i)
      if (origContractedPinList[i] == -1)
        origContractedPinList[i] = nonLocalMatches[localPins[i]];
  }

  /* send coarse hyperedges to appropriate processors via hash function */

  for (i = 0; i < numLocalHedges; ++i) {
    j = hEdgeOffsets[i + 1] - hEdgeOffsets[i];
    Funct::qsort(0, j - 1, &origContractedPinList[hEdgeOffsets[i]]);
  }

  contractedPinListLen = 0;
  numContractedHedges = 0;

  for (i = 0; i < numLocalHedges; ++i) {
    contractedHedgeOffsets.assign(numContractedHedges, contractedPinListLen);
    endOffset = hEdgeOffsets[i + 1];
    startOffset = hEdgeOffsets[i];

    contractedPinList.assign(contractedPinListLen,
                             origContractedPinList[startOffset]);
    contractedHedgeLen = 1;

    for (j = startOffset + 1; j < endOffset; ++j) {
      if (origContractedPinList[j] !=
          contractedPinList[contractedPinListLen + contractedHedgeLen - 1]) {
        contractedPinList.assign(contractedPinListLen + (contractedHedgeLen++),
                                 origContractedPinList[j]);
      }
    }

    if (contractedHedgeLen > 1) {
      contractedPinListLen += contractedHedgeLen;
      contractedHedgeWts.assign(numContractedHedges++, hEdgeWeights[i]);

      if (contractedHedgeLen > maxLocalCoarseHedgeLen)
        maxLocalCoarseHedgeLen = contractedHedgeLen;
    }
  }

  contractedHedgeOffsets.assign(numContractedHedges, contractedPinListLen);

  // ###
  // - compute maxCoarseHedgeLen and set in Funct
  // - compute hash-keys for each hyperedge
  // - send hyperedges to processor determined by
  //   corresponding hash key
  // ###

  MPI_Allreduce(&maxLocalCoarseHedgeLen, &maxCoarseHedgeLen, 1, MPI_INT,
                MPI_MAX, comm);

  Funct::setMaxHedgeLen(maxCoarseHedgeLen);

  for (i = 0; i < processors_; ++i)
    send_lens_[i] = 0;

  for (i = 0; i < numContractedHedges; ++i) {
    startOffset = contractedHedgeOffsets[i];
    endOffset = contractedHedgeOffsets[i + 1];
    contractedHedgeLen = endOffset - startOffset;

    p = Mod(
        Funct::computeHash(&contractedPinList[startOffset], contractedHedgeLen),
        processors_);

    data_out_sets_[p]->assign(send_lens_[p]++, contractedHedgeLen + 2);
    data_out_sets_[p]->assign(send_lens_[p]++, contractedHedgeWts[i]);

    for (j = startOffset; j < endOffset; ++j) {
      data_out_sets_[p]->assign(send_lens_[p]++, contractedPinList[j]);
    }
  }

  // ###
  // compute number of elements to send to other procs
  // ###

  j = 0;
  for (i = 0; i < processors_; ++i) {
    send_displs_[i] = j;
    j += send_lens_[i];
  }

  send_array_.reserve(j);
  totalToSend = j;

  j = 0;
  for (p = 0; p < processors_; ++p) {
    array = data_out_sets_[p]->data();
    arrayLen = send_lens_[p];

    for (i = 0; i < arrayLen; ++i) {
      send_array_[j++] = array[i];
    }
  }

  MPI_Alltoall(send_lens_.data(), 1, MPI_INT, receive_lens_.data(), 1, MPI_INT,
               comm);

  // ###
  // compute number of elements to receive from other procs
  // ###

  j = 0;
  for (i = 0; i < processors_; ++i) {
    receive_displs_[i] = j;
    j += receive_lens_[i];
  }

  receive_array_.reserve(j);
  totalToRecv = j;

  MPI_Alltoallv(send_array_.data(), send_lens_.data(),
                send_displs_.data(), MPI_INT, receive_array_.data(),
                receive_lens_.data(), receive_displs_.data(), MPI_INT, comm);

  // ###
  // now should have received all hyperedges destined for processor
  // now build the coarse hypergraph pin-list
  // ###

  dynamic_array<int> coarseLocalPins;
  dynamic_array<int> coarseHedgeOffsets;
  dynamic_array<int> coarseHedgeWts;

  int numCoarsePins = 0;
  int numCoarseHedges = 0;
  int coarseHedgeLen;
  int numSeen;
  int duplHedge;
  int length;
  int tryHedge;

  HashKey hashKey;

  ds::new_hyperedge_index_table table((int)(ceil((double)numLocalHedges * 1.5)));

  // ###
  // run through the recv data_
  // for each encountered hyperedge:
  // - compute hash-key
  // - check for duplicates in the hash table
  // - if duplicate not found, add to list
  // - if duplicate found, increment its weight
  // ###

  i = 0;
  coarseHedgeOffsets.assign(numCoarseHedges, numCoarsePins);

  // write_log(rank_, "entering hyperedges into table, totalToRecv = %d",
  // totalToRecv);

  while (i < totalToRecv) {
    /// if(rank_ == 10 && i > 645915)
    // write_log(rank_, "i = %d", i);

    coarseHedgeLen = receive_array_[i] - 2;
    hashKey = Funct::computeHash(&receive_array_[i + 2], coarseHedgeLen);

    /* find duplicate hyperedge (if exists) */

    numSeen = -1;
    duplHedge = -1;

#ifdef DEBUG_HYPERGRAPH
    int numSearches = 0;
#endif

    do {
      // if(rank_ == 10 && i > 645915)
      // write_log(rank_, "i = %d, [loop]", i);

      tryHedge = table.getHedgeIndex(hashKey, numSeen);

      // if(rank_ == 10 && i > 645915)
      // write_log(rank_, "i = %d, tryHedge = %d", i, tryHedge);

      if (tryHedge >= 0) {
#ifdef DEBUG_HYPERGRAPH
        assert(tryHedge < numCoarseHedges);
#endif
        startOffset = coarseHedgeOffsets[tryHedge];
        endOffset = coarseHedgeOffsets[tryHedge + 1];
        length = endOffset - startOffset;
#ifdef DEBUG_HYPERGRAPH
        assert(capacity_ > 0 && capacity_ < 0xFFFFFFF);
#endif
        if (length == coarseHedgeLen) {
          for (j = 0; j < length; ++j) {
#ifdef DEBUG_HYPERGRAPH
            assert(i + 2 + j < totalToRecv);
            assert(startOffset + j < numCoarsePins);
#endif
            if (receive_array_[i + 2 + j] != coarseLocalPins[startOffset + j])
              break;
          }
          if (j == length) {
            duplHedge = tryHedge;
            break;
          }
        }
      }

#ifdef DEBUG_HYPERGRAPH
      numSearches++;
      // if((numSearches % 100 == 0) && rank_ == 10)
      // write_log(rank_, "numSearches = %d", numSearches);
      assert(numSearches < 0xFFF);
#endif
    } while (numSeen >= 0);

    // if(rank_ == 10 && i > 645915)
    // write_log(rank_, "i = %d, duplHedge = %d", i, duplHedge);

    if (duplHedge == -1) {
      table.insertKey(hashKey, numCoarseHedges);
      coarseHedgeWts.assign(numCoarseHedges++, receive_array_[i + 1]);

      for (j = 0; j < coarseHedgeLen; ++j)
        coarseLocalPins.assign(numCoarsePins + j, receive_array_[i + 2 + j]);

      numCoarsePins += coarseHedgeLen;
      coarseHedgeOffsets.assign(numCoarseHedges, numCoarsePins);
    } else {
#ifdef DEBUG_HYPERGRAPH
      assert(duplHedge >= 0);
      assert(duplHedge < numCoarseHedges);
      assert(i + 1 < totalToRecv);
#endif
      coarseHedgeWts[duplHedge] += receive_array_[i + 1];
    }
#ifdef DEBUG_HYPERGRAPH
    assert(receive_array_[i] > 0);
#endif

    // if(i > 645915 && rank_ == 10)
    // write_log(rank_, "i = %d, receive_array_[i] = %d",i,receive_array_[i]);

    i += receive_array_[i];
  }

// write_log(rank_, "entered hyperedges into table");

#ifdef MEM_OPT
  coarseHedgeOffsets.reserve(numCoarseHedges + 1);
  coarseLocalPins.reserve(numCoarsePins);
  coarseHedgeWts.reserve(numCoarseHedges);
#endif

  /* now set the coarse hypergraph */

  coarse.setNumLocalHedges(numCoarseHedges);
  coarse.setNumLocalPins(numCoarsePins);
  coarse.allocHedgeMem(numCoarseHedges, numCoarsePins);

  int *cHedgeOffsets = coarse.getHedgeOffsetsArray();
  int *cPins = coarse.getLocalPinsArray();
  int *cHedgeWeights = coarse.getHedgeWeightsArray();

  for (i = 0; i < numCoarseHedges; ++i) {
    cHedgeWeights[i] = coarseHedgeWts[i];
    cHedgeOffsets[i] = coarseHedgeOffsets[i];
  }

  cHedgeOffsets[i] = coarseHedgeOffsets[i];

  for (i = 0; i < numCoarsePins; ++i)
    cPins[i] = coarseLocalPins[i];
}

void ParaHypergraph::contractRestrHyperedges(ParaHypergraph &coarse,
                                             MPI_Comm comm) {
  int maxLocalVertex = minVertexIndex + numLocalVertices;
  int totalToSend;
  int totalToRecv;
  int arrayLen;
  int p;
  int vertex;
  int startOffset;
  int endOffset;
  int contractedPinListLen;
  int contractedHedgeLen;
  int numContractedHedges;
  int maxLocalCoarseHedgeLen = 0;
  int maxCoarseHedgeLen;

  int *array;

  int i;
  int j;

#ifdef DEBUG_HYPERGRAPH
  int numTotCoarseVerts = coarse.getNumTotalVertices();
  for (i = 0; i < numLocalVertices; ++i)
    assert(matchVector[i] >= 0 && matchVector[i] < numTotCoarseVerts);
#endif

  dynamic_array<int> minFineIdxOnProc(processors_);
  dynamic_array<int> procs(processors_);
  dynamic_array<int> contractedPinList;
  dynamic_array<int> contractedHedgeOffsets;
  dynamic_array<int> contractedHedgeWts;
  dynamic_array<int> origContractedPinList(numLocalPins);
  dynamic_array<int> copyOfReq;

  MPI_Allgather(&minVertexIndex, 1, MPI_INT, minFineIdxOnProc.data(), 1,
                MPI_INT, comm);

  for (i = 0; i < processors_; ++i)
    procs[i] = i;

  ds::complete_binary_tree <int> vFineToProc(procs.data(),
                                             minFineIdxOnProc.data(), processors_);

  for (i = 0; i < processors_; ++i)
    send_lens_[i] = 0;

  ds::bit_field sentRequests(numTotalVertices);
  sentRequests.unset();

  // ###
  // compute all the requests for remote vertex matches
  // ###

  for (i = 0; i < numLocalPins; ++i) {
    vertex = localPins[i];

    if (vertex < minVertexIndex || vertex >= maxLocalVertex) {
      if (!sentRequests(vertex)) {
        p = vFineToProc.root_value(vertex);
        data_out_sets_[p]->assign(send_lens_[p]++, vertex);
        sentRequests.set(vertex);
      }

      origContractedPinList[i] = -1;
    } else {
      origContractedPinList[i] = matchVector[vertex - minVertexIndex];
    }
  }

  // ###
  // compute number of elements to send to other procs
  // ###

  j = 0;
  for (i = 0; i < processors_; ++i) {
    send_displs_[i] = j;
    j += send_lens_[i];
  }

  send_array_.reserve(j);
  copyOfReq.reserve(j);
  totalToSend = j;

  j = 0;
  for (p = 0; p < processors_; ++p) {
    array = data_out_sets_[p]->data();
    arrayLen = send_lens_[p];

    for (i = 0; i < arrayLen; ++i) {
      vertex = array[i];
      send_array_[j] = vertex;
      copyOfReq[j++] = vertex;
    }
  }

  MPI_Alltoall(send_lens_.data(), 1, MPI_INT, receive_lens_.data(), 1, MPI_INT,
               comm);

  // ###
  // compute number of elements to receive from other procs
  // ###

  j = 0;
  for (i = 0; i < processors_; ++i) {
    receive_displs_[i] = j;
    j += receive_lens_[i];
  }

  receive_array_.reserve(j);
  totalToRecv = j;

  MPI_Alltoallv(send_array_.data(), send_lens_.data(),
                send_displs_.data(), MPI_INT, receive_array_.data(),
                receive_lens_.data(), receive_displs_.data(), MPI_INT, comm);

  // ###
  // now have received all requests and sent out our requests
  // the reply communication will have the dual dimensions
  // ###

  send_array_.reserve(totalToRecv);

  for (i = 0; i < totalToRecv; ++i) {
#ifdef DEBUG_HYPERGRAPH
    assert(receive_array_[i] >= minVertexIndex &&
           receive_array_[i] < maxLocalVertex);
#endif
    send_array_[i] = matchVector[receive_array_[i] - minVertexIndex];
  }

  receive_array_.reserve(totalToSend);

  MPI_Alltoallv(send_array_.data(), receive_lens_.data(),
                receive_displs_.data(), MPI_INT, receive_array_.data(),
                send_lens_.data(), send_displs_.data(), MPI_INT, comm);

  // ###
  // now the requested vertices are in the copyOfReq data_
  // while their corresponding matchVector values are in
  // the corresponding location in the recvArray
  // ###

  if (numLocalPins < numTotalVertices / 2) {
    ds::map_from_pos_int<int> storedRequests(numLocalPins);

    for (i = 0; i < totalToSend; ++i)
      storedRequests.insert(copyOfReq[i], receive_array_[i]);

    /* contract remaining local pins */

    for (i = 0; i < numLocalPins; ++i)
      if (origContractedPinList[i] == -1)
        origContractedPinList[i] = storedRequests.get(localPins[i]);
  } else {
    dynamic_array<int> nonLocalMatches(numTotalVertices);

    for (i = 0; i < totalToSend; ++i)
      nonLocalMatches[copyOfReq[i]] = receive_array_[i];

    /* contract remaining local pins */

    for (i = 0; i < numLocalPins; ++i)
      if (origContractedPinList[i] == -1)
        origContractedPinList[i] = nonLocalMatches[localPins[i]];
  }

  /* send coarse hyperedges to appropriate processors via hash function */

  for (i = 0; i < numLocalHedges; ++i) {
    j = hEdgeOffsets[i + 1] - hEdgeOffsets[i];
    Funct::qsort(0, j - 1, &origContractedPinList[hEdgeOffsets[i]]);
  }

  contractedPinListLen = 0;
  numContractedHedges = 0;

  for (i = 0; i < numLocalHedges; ++i) {
    contractedHedgeOffsets.assign(numContractedHedges, contractedPinListLen);
    endOffset = hEdgeOffsets[i + 1];
    startOffset = hEdgeOffsets[i];

    contractedPinList.assign(contractedPinListLen,
                             origContractedPinList[startOffset]);
    contractedHedgeLen = 1;

    for (j = startOffset + 1; j < endOffset; ++j) {
      if (origContractedPinList[j] !=
          contractedPinList[contractedPinListLen + contractedHedgeLen - 1]) {
        contractedPinList.assign(contractedPinListLen + (contractedHedgeLen++),
                                 origContractedPinList[j]);
      }
    }

    if (contractedHedgeLen > 1) {
      contractedPinListLen += contractedHedgeLen;
      contractedHedgeWts.assign(numContractedHedges++, hEdgeWeights[i]);

      if (contractedHedgeLen > maxLocalCoarseHedgeLen)
        maxLocalCoarseHedgeLen = contractedHedgeLen;
    }
  }

  contractedHedgeOffsets.assign(numContractedHedges, contractedPinListLen);

  // ###
  // - compute maxCoarseHedgeLen and set in Funct
  // - compute hash-keys for each hyperedge
  // - send hyperedges to processor determined by
  //   corresponding hash key
  // ###

  MPI_Allreduce(&maxLocalCoarseHedgeLen, &maxCoarseHedgeLen, 1, MPI_INT,
                MPI_MAX, comm);

  Funct::setMaxHedgeLen(maxCoarseHedgeLen);

  for (i = 0; i < processors_; ++i)
    send_lens_[i] = 0;

  for (i = 0; i < numContractedHedges; ++i) {
    startOffset = contractedHedgeOffsets[i];
    endOffset = contractedHedgeOffsets[i + 1];
    contractedHedgeLen = endOffset - startOffset;

    p = Mod(
        Funct::computeHash(&contractedPinList[startOffset], contractedHedgeLen),
        processors_);

    data_out_sets_[p]->assign(send_lens_[p]++, contractedHedgeLen + 2);
    data_out_sets_[p]->assign(send_lens_[p]++, contractedHedgeWts[i]);

    for (j = startOffset; j < endOffset; ++j) {
      data_out_sets_[p]->assign(send_lens_[p]++, contractedPinList[j]);
    }
  }

  // ###
  // compute number of elements to send to other procs
  // ###

  j = 0;
  for (i = 0; i < processors_; ++i) {
    send_displs_[i] = j;
    j += send_lens_[i];
  }

  send_array_.reserve(j);
  totalToSend = j;

  j = 0;
  for (p = 0; p < processors_; ++p) {
    array = data_out_sets_[p]->data();
    arrayLen = send_lens_[p];

    for (i = 0; i < arrayLen; ++i) {
      send_array_[j++] = array[i];
    }
  }

  MPI_Alltoall(send_lens_.data(), 1, MPI_INT, receive_lens_.data(), 1, MPI_INT,
               comm);

  // ###
  // compute number of elements to receive from other procs
  // ###

  j = 0;
  for (i = 0; i < processors_; ++i) {
    receive_displs_[i] = j;
    j += receive_lens_[i];
  }

  receive_array_.reserve(j);
  totalToRecv = j;

  MPI_Alltoallv(send_array_.data(), send_lens_.data(),
                send_displs_.data(), MPI_INT, receive_array_.data(),
                receive_lens_.data(), receive_displs_.data(), MPI_INT, comm);

  // ###
  // now should have received all hyperedges destined for processor
  // now build the coarse hypergraph pin-list
  // ###

  dynamic_array<int> coarseLocalPins;    // = new dynamic_array<int>(2048);
  dynamic_array<int> coarseHedgeOffsets; // = new dynamic_array<int>(1024);
  dynamic_array<int> coarseHedgeWts;     // = new dynamic_array<int>(1024);

  int numCoarsePins = 0;
  int numCoarseHedges = 0;
  int coarseHedgeLen;
  int numSeen;
  int duplHedge;
  int length;
  int tryHedge;

  HashKey hashKey;

  ds::new_hyperedge_index_table table((int)(ceil((double)numLocalHedges * 1.1)));

  // ###
  // run through the recv data_
  // for each encountered hyperedge:
  // - compute hash-key
  // - check for duplicates in the hash table
  // - if duplicate not found, add to list
  // - if duplicate found, increment its weight
  // ###

  /* NEED TO MODIFY THE RESTR HEDGE CONTRACTION TO CORRESPOND WITH NORMAL */

  i = 0;
  coarseHedgeOffsets.assign(numCoarseHedges, numCoarsePins);

  while (i < totalToRecv) {
    coarseHedgeLen = receive_array_[i] - 2;
    hashKey = Funct::computeHash(&receive_array_[i + 2], coarseHedgeLen);

    // ###
    // find duplicate hyperedge (if exists)
    // ###

    numSeen = -1;
    duplHedge = -1;

    do {
      tryHedge = table.getHedgeIndex(hashKey, numSeen);

      if (tryHedge >= 0) {
#ifdef DEBUG_HYPERGRAPH
        assert(tryHedge < numCoarseHedges);
#endif
        startOffset = coarseHedgeOffsets[tryHedge];
        endOffset = coarseHedgeOffsets[tryHedge + 1];
        length = endOffset - startOffset;

        if (length == coarseHedgeLen) {
          for (j = 0; j < length; ++j) {
#ifdef DEBUG_HYPERGRAPH
            assert(i + 2 + j < totalToRecv);
            assert(startOffset + j < numCoarsePins);
#endif
            if (receive_array_[i + 2 + j] != coarseLocalPins[startOffset + j])
              break;
          }
          if (j == length) {
            duplHedge = tryHedge;
            break;
          }
        }
      }
    } while (numSeen >= 0);

    if (duplHedge == -1) {
      table.insertKey(hashKey, numCoarseHedges);
      coarseHedgeWts.assign(numCoarseHedges++, receive_array_[i + 1]);

      for (j = 0; j < coarseHedgeLen; ++j)
        coarseLocalPins.assign(numCoarsePins + j, receive_array_[i + 2 + j]);

      numCoarsePins += coarseHedgeLen;
      coarseHedgeOffsets.assign(numCoarseHedges, numCoarsePins);
    } else {
#ifdef DEBUG_HYPERGRAPH
      assert(duplHedge >= 0);
      assert(duplHedge < numCoarseHedges);
      assert(i + 1 < totalToRecv);
#endif
      coarseHedgeWts[duplHedge] += receive_array_[i + 1];
    }

    i += receive_array_[i];
  }

#ifdef MEM_OPT
  coarseHedgeOffsets.reserve(numCoarseHedges + 1);
  coarseLocalPins.reserve(numCoarsePins);
  coarseHedgeWts.reserve(numCoarseHedges);
#endif

  // ###
  // now set the coarse hypergraph
  // ###

  coarse.setNumLocalHedges(numCoarseHedges);
  coarse.setNumLocalPins(numCoarsePins);
  coarse.allocHedgeMem(numCoarseHedges, numCoarsePins);

  int *cHedgeOffsets = coarse.getHedgeOffsetsArray();
  int *cPins = coarse.getLocalPinsArray();
  int *cHedgeWeights = coarse.getHedgeWeightsArray();

  for (i = 0; i < numCoarseHedges; ++i) {
    cHedgeWeights[i] = coarseHedgeWts[i];
    cHedgeOffsets[i] = coarseHedgeOffsets[i];
  }

  cHedgeOffsets[i] = coarseHedgeOffsets[i];

  for (i = 0; i < numCoarsePins; ++i)
    cPins[i] = coarseLocalPins[i];
}

void ParaHypergraph::projectPartitions(ParaHypergraph &coarse, MPI_Comm comm) {
  int totCoarseV = coarse.getNumTotalVertices();
  int numLocCoarseV = coarse.getNumLocalVertices();
  int minCoarseVindex = coarse.getMinVertexIndex();
  int maxCoarseVindex = minCoarseVindex + numLocCoarseV;
  int numCoarsePartitions = coarse.getNumPartitions();

  int *coarsePartitionVector = coarse.getPartitionArray();
  int *coarsePartitionOffsets = coarse.getPartitionOffsetsArray();
  int *coarsePartitionCutsizes = coarse.getCutsizesArray();
  int *array;

  int vCoarsePerProc = totCoarseV / processors_;
  int vertex;
  int proc;
  int totToSend;
  int totToRecv;
  int sendLength;
  int numRequestingLocalVerts;
  int locCoarseVindex;

  int i;
  int j;
  int ij;

  dynamic_array<int> requestingLocalVerts;

  for (i = 0; i < processors_; ++i)
    send_lens_[i] = 0;

  // ###
  // initialise the local partition structures
  // ###

  numPartitions = numCoarsePartitions;

  partitionCutsizesVector.reserve(numPartitions);
  partitionOffsetsVector.reserve(numPartitions + 1);
  partitionVector.reserve(numPartitions * numLocalVertices);

  for (i = 0; i < numPartitions; ++i) {
#ifdef DEBUG_HYPERGRAPH
    assert(coarsePartitionCutsizes[i] > 0);
#endif
    partitionCutsizesVector[i] = coarsePartitionCutsizes[i];
  }

  j = 0;

  for (i = 0; i <= numPartitions; ++i) {
    partitionOffsetsVector[i] = j;
    j += numLocalVertices;
  }

#ifdef DEBUG_HYPERGRAPH
  for (i = 0; i < partitionOffsetsVector[numPartitions]; ++i)
    partitionVector[i] = -1;
#endif

  for (i = 0; i < numLocalVertices; ++i) {
#ifdef DEBUG_HYPERGRAPH
    assert(matchVector[i] >= 0 && matchVector[i] < totCoarseV);
#endif

    ij = matchVector[i];

    if (ij >= minCoarseVindex && ij < maxCoarseVindex) {
      locCoarseVindex = ij - minCoarseVindex;

      for (j = 0; j < numPartitions; ++j) {
        partitionVector[partitionOffsetsVector[j] + i] =
            coarsePartitionVector[coarsePartitionOffsets[j] + locCoarseVindex];
      }
    } else {
      proc = std::min(ij / vCoarsePerProc, processors_ - 1);

      data_out_sets_[proc]->assign(send_lens_[proc]++, i);
      data_out_sets_[proc]->assign(send_lens_[proc]++, ij);
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

#ifdef DEBUG_HYPERGRAPH
  totToSend = ij;
#endif

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

#ifdef DEBUG_HYPERGRAPH
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
  // local vertex partitions
  // dimensions of cummunication
  // are reversed!
  // ###

  totToSend = totToRecv * numPartitions;
  send_array_.reserve(totToSend);

  ij = 0;

  for (i = 0; i < totToRecv; ++i) {
#ifdef DEBUG_HYPERGRAPH
    assert(receive_array_[i] >= minCoarseVindex &&
           receive_array_[i] < maxCoarseVindex);
#endif

    vertex = receive_array_[i] - minCoarseVindex;

    for (j = 0; j < numPartitions; ++j) {
      send_array_[ij++] =
          coarsePartitionVector[coarsePartitionOffsets[j] + vertex];
    }
  }

  for (i = 0; i < processors_; ++i) {
    ij = receive_lens_[i];
    receive_lens_[i] = send_lens_[i] * numPartitions;
    send_lens_[i] = ij * numPartitions;
  }

  ij = 0;
  j = 0;

  for (i = 0; i < processors_; ++i) {
    send_displs_[i] = ij;
    receive_displs_[i] = j;

    ij += send_lens_[i];
    j += receive_lens_[i];
  }

  receive_array_.reserve(j);

#ifdef DEBUG_HYPERGRAPH
  assert(send_array_.getLength() == ij);
  assert(numRequestingLocalVerts * numPartitions == j);
  assert(requestingLocalVerts.getLength() * numPartitions == j);
#endif

  MPI_Alltoallv(send_array_.data(), send_lens_.data(),
                send_displs_.data(), MPI_INT, receive_array_.data(),
                receive_lens_.data(), receive_displs_.data(), MPI_INT, comm);

  // ###
  // finish off initialising the
  // partition vector
  // ###

  ij = 0;

  for (i = 0; i < numRequestingLocalVerts; ++i) {
    vertex = requestingLocalVerts[i];

    for (j = 0; j < numPartitions; ++j) {
      partitionVector[partitionOffsetsVector[j] + vertex] = receive_array_[ij++];
    }
  }

#ifdef DEBUG_HYPERGRAPH
  for (i = 0; i < partitionOffsetsVector[numPartitions]; ++i)
    assert(partitionVector[i] != -1);
#endif
}

void ParaHypergraph::resetVectors() {
  int i;

  for (i = 0; i < numLocalVertices; ++i)
    matchVector[i] = -1;

  numPartitions = 0;

  vToOrigV.reserve(0);
  partitionVector.reserve(0);
  partitionOffsetsVector.reserve(0);
  partitionCutsizesVector.reserve(0);

  free_memory();
}

void ParaHypergraph::removeBadPartitions(double cutThreshold) {
#ifdef DEBUG_HYPERGRAPH
  assert(numPartitions > 0);
#endif

  int i;
  int j;

  int bestPartition = 0;
  int bestCut = partitionCutsizesVector[0];
  int acceptedCut;
  int diffInCut;
  int indexIntoOld;
  int indexIntoNew;
  int pSeenBefore;
  int endOffset;
  int numNewPartitions = 0;

  for (i = 1; i < numPartitions; ++i) {
    if (partitionCutsizesVector[i] < bestCut) {
      bestCut = partitionCutsizesVector[i];
      bestPartition = i;
    }
  }

  diffInCut =
      static_cast<int>(floor(static_cast<float>(bestCut) * cutThreshold));
  acceptedCut = bestCut + diffInCut;

  indexIntoNew = 0;
  indexIntoOld = 0;

  for (i = 0; i < numPartitions; ++i) {
    if (partitionCutsizesVector[i] <= acceptedCut) {
      pSeenBefore = 0;

      for (j = 0; j < numNewPartitions; ++j) {
        if (partitionCutsizesVector[j] == partitionCutsizesVector[i])
          pSeenBefore = 1;
      }
      if (pSeenBefore == 0) {
        if (indexIntoOld > indexIntoNew) {
          endOffset = partitionOffsetsVector[i + 1];

          for (; indexIntoOld < endOffset; ++indexIntoOld) {
            partitionVector[indexIntoNew++] = partitionVector[indexIntoOld];
          }
        } else {
          indexIntoOld += numLocalVertices;
          indexIntoNew += numLocalVertices;
        }

        partitionCutsizesVector[numNewPartitions++] =
            partitionCutsizesVector[i];
      } else {
        indexIntoOld += numLocalVertices;
      }
    } else {
      indexIntoOld += numLocalVertices;
    }
  }

  numPartitions = numNewPartitions;
}

void ParaHypergraph::setNumberPartitions(int nP) {
  int i;
  int j = 0;

  numPartitions = nP;

  partitionCutsizesVector.reserve(nP);
  partitionOffsetsVector.reserve(nP + 1);

  for (i = 0; i <= numPartitions; ++i) {
    partitionOffsetsVector[i] = j;
    j += numLocalVertices;
  }

  partitionVector.reserve(partitionOffsetsVector[numPartitions]);
}

void ParaHypergraph::computePartitionChars(int pNum, int numParts,
                                           double constraint, std::ostream &out,
                                           MPI_Comm comm) {
#ifdef DEBUG_HYPERGRAPH
  assert(numPartitions > 0);
  assert(pNum < numPartitions);
#endif

  int i;

  int cut;
  int maxAllowedPartWt;
  int maxPartWt;
  int minPartWt;
  int totHypergraphWeight;

  int *pVector = &partitionVector[partitionOffsetsVector[pNum]];

  double avePartWt;

  dynamic_array<int> locPartWeights(numParts);
  dynamic_array<int> partWeights(numParts);

  cut = calcCutsize(numParts, pNum, comm);

  for (i = 0; i < numParts; ++i)
    locPartWeights[i] = 0;

  for (i = 0; i < numLocalVertices; ++i)
    locPartWeights[pVector[i]] += vWeight[i];

  MPI_Reduce(&localVertexWt, &totHypergraphWeight, 1, MPI_INT, MPI_SUM, 0,
             comm);
  MPI_Reduce(locPartWeights.data(), partWeights.data(), numParts,
             MPI_INT, MPI_SUM, 0, comm);

  if (rank_ == 0) {
    avePartWt = static_cast<double>(totHypergraphWeight) / numParts;
    maxAllowedPartWt =
        static_cast<int>(floor(avePartWt + avePartWt * constraint));

    maxPartWt = 0;
    minPartWt = LARGE_CONSTANT;

    for (i = 0; i < numParts; ++i) {
      if (partWeights[i] > maxPartWt)
        maxPartWt = partWeights[i];
      if (partWeights[i] < minPartWt)
        minPartWt = partWeights[i];
    }

    out << "****** partition summary ******" << std::endl
        << std::endl
        << "\tcut = " << cut << std::endl
        << "\tmaxAllowedPartWt = " << maxAllowedPartWt << std::endl
        << "\tminPartWt = " << minPartWt << std::endl
        << "\tmaxPartWt = " << maxPartWt << std::endl;
  }

  MPI_Barrier(comm);
}

void ParaHypergraph::copyInPartition(const int *partition, int numV,
                                     int nP) {
#ifdef DEBUG_HYPERGRAPH
  assert(numLocalVertices == numV);
  assert(nP >= 0 && nP < numPartitions);
#endif

  int i;

  int endOffset = partitionOffsetsVector[nP + 1];
  int startOffset = partitionOffsetsVector[nP];

  for (i = startOffset; i < endOffset; ++i)
    partitionVector[i] = partition[i - startOffset];
}

void ParaHypergraph::copyOutPartition(int *partition, int numV,
                                      int nP) const {
#ifdef DEBUG_HYPERGRAPH
  assert(numLocalVertices == numV);
  assert(nP >= 0 && nP < numPartitions);
#endif

  int i;

  int startOffset = partitionOffsetsVector[nP];
  int endOffset = partitionOffsetsVector[nP + 1];

  for (i = startOffset; i < endOffset; ++i)
    partition[i] = partitionVector[i - startOffset];
}

int ParaHypergraph::keepBestPartition() {
  int i;
  int bestOffset;

  int bestPartition = 0;
  int bestCut = partitionCutsizesVector[0];

  for (i = 1; i < numPartitions; ++i) {
    if (partitionCutsizesVector[i] < bestCut) {
      bestPartition = i;
      bestCut = partitionCutsizesVector[i];
    }
  }

  if (bestPartition != 0) {
    bestOffset = partitionOffsetsVector[bestPartition];

    for (i = 0; i < numLocalVertices; ++i) {
      partitionVector[i] = partitionVector[bestOffset + i];
    }
  }

  partitionCutsizesVector[0] = bestCut;
  numPartitions = 1;

  return bestCut;
}

void ParaHypergraph::prescribedVertexShuffle(int *mapToOrigV, int *prescArray,
                                             MPI_Comm comm) {
  prescribedVertexShuffle(prescArray, numLocalVertices, comm);
  shiftVerticesToBalance(comm);

  int vPerProc = numTotalVertices / processors_;
  int maxLocalVertex = minVertexIndex + numLocalVertices;
  int totalToRecv;
  int totalToSend;
  int arrayLen;
  int vertex;
  int *array1;
  int *array2;

  int i;
  int j;
  int ij;

  dynamic_array<int> copyOfMapToOrigV(numLocalVertices);
  dynamic_array<int> askingVertex;

  dynamic_array<dynamic_array<int> *> askingVertices(processors_);

  for (i = 0; i < numLocalVertices; ++i)
    copyOfMapToOrigV[i] = mapToOrigV[i];

  for (i = 0; i < processors_; ++i) {
    send_lens_[i] = 0;
    askingVertices[i] = new dynamic_array<int>(0);
  }

  /* compute mapToOrigV entries required from un-shuffled hypergraph */

  for (i = 0; i < numLocalVertices; ++i) {
    vertex = vToOrigV[i];

    if (vertex < minVertexIndex || vertex >= maxLocalVertex) {
      j = std::min(vertex / vPerProc, processors_ - 1);

      askingVertices[j]->assign(send_lens_[j], i);
      data_out_sets_[j]->assign(send_lens_[j]++, vertex);
    } else {
      mapToOrigV[i] = copyOfMapToOrigV[vertex - minVertexIndex];
    }
  }

  /* compute number of elements to send to other procs */

  j = 0;
  for (i = 0; i < processors_; ++i) {
    send_displs_[i] = j;
    j += send_lens_[i];
  }

  send_array_.reserve(j);
  askingVertex.reserve(j);
  totalToSend = j;

  j = 0;
  for (ij = 0; ij < processors_; ++ij) {
    array1 = data_out_sets_[ij]->data();
    array2 = askingVertices[ij]->data();

    arrayLen = send_lens_[ij];

    for (i = 0; i < arrayLen; ++i) {
      send_array_[j] = array1[i];
      askingVertex[j++] = array2[i];
    }
  }

#ifdef DEBUG_HYPERGRAPH
  assert(j == totalToSend);
#endif

  MPI_Alltoall(send_lens_.data(), 1, MPI_INT, receive_lens_.data(), 1, MPI_INT,
               comm);

  /* compute number of elements to receive from other procs */

  j = 0;
  for (i = 0; i < processors_; ++i) {
    receive_displs_[i] = j;
    j += receive_lens_[i];
  }

  receive_array_.reserve(j);
  totalToRecv = j;

  MPI_Alltoallv(send_array_.data(), send_lens_.data(),
                send_displs_.data(), MPI_INT, receive_array_.data(),
                receive_lens_.data(), receive_displs_.data(), MPI_INT, comm);

  /*
    now have received all requests and sent out our requests
    the reply communication will have the dual dimensions
  */

  send_array_.reserve(totalToRecv);

  for (i = 0; i < totalToRecv; ++i) {
#ifdef DEBUG_HYPERGRAPH
    assert(receive_array_[i] >= minVertexIndex &&
           receive_array_[i] < maxLocalVertex);
#endif
    send_array_[i] = copyOfMapToOrigV[receive_array_[i] - minVertexIndex];
#ifdef DEBUG_HYPERGRAPH
    assert(send_array_[i] >= 0 && send_array_[i] < numTotalVertices);
#endif
  }

  receive_array_.reserve(totalToSend);

  MPI_Alltoallv(send_array_.data(), receive_lens_.data(),
                receive_displs_.data(), MPI_INT, receive_array_.data(),
                send_lens_.data(), send_displs_.data(), MPI_INT, comm);

  /*
     now the requested vertices are in the copyOfReq data_
     while their corresponding matchVector values are in
     the corresponding location in the receive_array_
  */

  for (i = 0; i < totalToSend; ++i) {
#ifdef DEBUG_HYPERGRAPH
    assert(receive_array_[i] >= 0 && receive_array_[i] < numTotalVertices);
#endif
    mapToOrigV[askingVertex[i]] = receive_array_[i];
  }

  for (i = 0; i < processors_; ++i)
    DynaMem::deletePtr<dynamic_array<int> >(askingVertices[i]);
}

void ParaHypergraph::prescribedVertexShuffle(int *prescribedAssignment,
                                             int nLocVer, MPI_Comm comm) {
#ifdef DEBUG_HYPERGRAPH
  assert(numLocalVertices == nLocVer);
#endif

  int i;

  dynamic_array<int> localVPerProc(processors_);

  for (i = 0; i < processors_; ++i)
    localVPerProc[i] = 0;

  for (i = 0; i < numLocalVertices; ++i) {
#ifdef DEBUG_HYPERGRAPH
    assert(prescribedAssignment[i] >= 0 && prescribedAssignment[i] < processors_);
#endif
    ++localVPerProc[prescribedAssignment[i]];
  }

  shuffleVertices(prescribedAssignment, localVPerProc.data(), comm);
}

void ParaHypergraph::shuffleVerticesByPartition(int nParts, MPI_Comm comm) {
#ifdef DEBUG_HYPERGRAPH
  assert(numPartitions > 0);
  assert(Mod(nParts, processors_) == 0);
#endif

  int i;
  int j;

  int numPartsPerProc = nParts / processors_;

  dynamic_array<int> vToProc(numLocalVertices);
  dynamic_array<int> localVPerProc(processors_);

  for (i = 0; i < processors_; ++i)
    localVPerProc[i] = 0;

  for (i = 0; i < numLocalVertices; ++i) {
    j = partitionVector[i] / numPartsPerProc;
    vToProc[i] = j;
    ++localVPerProc[j];
  }

  shuffleVertices(vToProc.data(), localVPerProc.data(), comm);
}

void ParaHypergraph::randomVertexShuffle(MPI_Comm comm) {
  int numVerticesEvenlyAllocated;
  int numSpareVertices;
  int totSpareVertices;
  int toProc;

  int i;
  int j;
  int ij;

  dynamic_array<int> vertices(numLocalVertices);
  dynamic_array<int> vToProc(numLocalVertices);
  dynamic_array<int> localVPerProc(processors_);
  dynamic_array<int> vSpareToProc;
  dynamic_array<int> indexIntoSpares(processors_);

  for (i = 0; i < processors_; ++i)
    localVPerProc[i] = 0;

  for (i = 0; i < numLocalVertices; ++i)
    vertices[i] = i;

  Funct::randomPermutation(vertices.data(), numLocalVertices);

  numSpareVertices = numLocalVertices % processors_;
  numVerticesEvenlyAllocated = numLocalVertices - numSpareVertices;

  MPI_Allgather(&numSpareVertices, 1, MPI_INT, indexIntoSpares.data(), 1,
                MPI_INT, comm);

  totSpareVertices = 0;

  for (i = 0; i < processors_; ++i) {
    j = totSpareVertices;
    totSpareVertices += indexIntoSpares[i];
    indexIntoSpares[i] = j;
  }

  vSpareToProc.reserve(totSpareVertices);

  j = totSpareVertices / processors_;
  ij = Mod(totSpareVertices, processors_);

  if (j == 0) {
    for (i = 0; i < totSpareVertices; ++i) {
      vSpareToProc[i] = processors_ - 1;
    }
  } else {
    for (i = 0; i < totSpareVertices - ij; ++i) {
      vSpareToProc[i] = Mod(i, processors_);
    }

    for (i = totSpareVertices - ij; i < totSpareVertices; ++i) {
      vSpareToProc[i] = processors_ - 1;
    }
  }

  for (i = 0; i < numVerticesEvenlyAllocated; ++i) {
    j = Mod(i, processors_);
    vToProc[vertices[i]] = j;
    ++localVPerProc[j];
  }

  for (; i < numLocalVertices; ++i) {
    j = i - numVerticesEvenlyAllocated;
    toProc = vSpareToProc[indexIntoSpares[rank_] + j];
    vToProc[vertices[i]] = toProc;
    ++localVPerProc[toProc];
  }

  shuffleVertices(vToProc.data(), localVPerProc.data(), comm);
}

void ParaHypergraph::randomVertexShuffle(int *mapToOrigV, MPI_Comm comm) {
  int numVerticesEvenlyAllocated;
  int numSpareVertices;
  int totSpareVertices;
  int toProc;

  int i;
  int j;
  int ij;

  dynamic_array<int> vertices(numLocalVertices);
  dynamic_array<int> vToProc(numLocalVertices);
  dynamic_array<int> localVPerProc(processors_);
  dynamic_array<int> vSpareToProc;
  dynamic_array<int> indexIntoSpares(processors_);

  /* first compute the V->proc map */

  for (i = 0; i < numLocalVertices; ++i)
    vertices[i] = i;

  Funct::randomPermutation(vertices.data(), numLocalVertices);

  numSpareVertices = Mod(numLocalVertices, processors_);
  numVerticesEvenlyAllocated = numLocalVertices - numSpareVertices;

  MPI_Allgather(&numSpareVertices, 1, MPI_INT, indexIntoSpares.data(), 1,
                MPI_INT, comm);

  totSpareVertices = 0;

  for (i = 0; i < processors_; ++i) {
    j = totSpareVertices;
    totSpareVertices += indexIntoSpares[i];
    indexIntoSpares[i] = j;
  }

  vSpareToProc.reserve(totSpareVertices);

  j = totSpareVertices / processors_;
  ij = Mod(totSpareVertices, processors_);

  if (j == 0) {
    for (i = 0; i < totSpareVertices; ++i) {
      vSpareToProc[i] = processors_ - 1;
    }
  } else {
    for (i = 0; i < totSpareVertices - ij; ++i) {
      vSpareToProc[i] = Mod(i, processors_);
    }

    for (i = totSpareVertices - ij; i < totSpareVertices; ++i) {
      vSpareToProc[i] = processors_ - 1;
    }
  }

  for (i = 0; i < numVerticesEvenlyAllocated; ++i) {
    j = Mod(i, processors_);
    vToProc[vertices[i]] = j;
  }

  j = numLocalVertices / processors_;
  for (i = 0; i < processors_; ++i)
    localVPerProc[i] = j;

  for (i = numVerticesEvenlyAllocated; i < numLocalVertices; ++i) {
    j = i - numVerticesEvenlyAllocated;
    toProc = vSpareToProc[indexIntoSpares[rank_] + j];
    vToProc[vertices[i]] = toProc;
    ++localVPerProc[toProc];
  }

  /* map V->proc is now computed */

  shuffleVerticesAftRandom(vToProc.data(), localVPerProc.data(),
                           mapToOrigV, comm);
}

void ParaHypergraph::randomVertexShuffle(ParaHypergraph &fG, MPI_Comm comm) {
  int numVerticesEvenlyAllocated;
  int numSpareVertices;
  int totSpareVertices;
  int toProc;

  int i;
  int j;
  int ij;

  dynamic_array<int> vertices(numLocalVertices);
  dynamic_array<int> vToProc(numLocalVertices);
  dynamic_array<int> localVPerProc(processors_);
  dynamic_array<int> vSpareToProc;
  dynamic_array<int> indexIntoSpares(processors_);

  /* first compute the V->proc map */

  for (i = 0; i < numLocalVertices; ++i)
    vertices[i] = i;

  Funct::randomPermutation(vertices.data(), numLocalVertices);

  numSpareVertices = Mod(numLocalVertices, processors_);
  numVerticesEvenlyAllocated = numLocalVertices - numSpareVertices;

  MPI_Allgather(&numSpareVertices, 1, MPI_INT, indexIntoSpares.data(), 1,
                MPI_INT, comm);

  totSpareVertices = 0;

  for (i = 0; i < processors_; ++i) {
    j = totSpareVertices;
    totSpareVertices += indexIntoSpares[i];
    indexIntoSpares[i] = j;
  }

  vSpareToProc.reserve(totSpareVertices);

  j = totSpareVertices / processors_;
  ij = Mod(totSpareVertices, processors_);

  if (j == 0) {
    for (i = 0; i < totSpareVertices; ++i) {
      vSpareToProc[i] = processors_ - 1;
    }
  } else {
    for (i = 0; i < totSpareVertices - ij; ++i) {
      vSpareToProc[i] = Mod(i, processors_);
    }

    for (i = totSpareVertices - ij; i < totSpareVertices; ++i) {
      vSpareToProc[i] = processors_ - 1;
    }
  }

  for (i = 0; i < numVerticesEvenlyAllocated; ++i) {
    j = Mod(i, processors_);
    vToProc[vertices[i]] = j;
  }

  j = numLocalVertices / processors_;
  for (i = 0; i < processors_; ++i)
    localVPerProc[i] = j;

  for (i = numVerticesEvenlyAllocated; i < numLocalVertices; ++i) {
    j = i - numVerticesEvenlyAllocated;
    toProc = vSpareToProc[indexIntoSpares[rank_] + j];
    vToProc[vertices[i]] = toProc;
    ++localVPerProc[toProc];
  }

  /* map V->proc is now computed */

  shuffleVerticesAftRandom(vToProc.data(), localVPerProc.data(), fG,
                           comm);
}

void ParaHypergraph::shuffleVertices(int *vToProc, int *localVPerProc,
                                     MPI_Comm comm) {
#ifdef DEBUG_HYPERGRAPH
  assert(numPartitions == 0 || numPartitions == 1);
#endif

  int i;
  int j;
  int ij;

  int maxLocalVertex = minVertexIndex + numLocalVertices;
  int vFinePerProc = numTotalVertices / processors_;
  int totalToSend;
  int totalToRecv;
  int vertex;
  int startOffset;
  int arrayLen;
  int *array;

  dynamic_array<int> oldIndexToNew(numLocalVertices);
  dynamic_array<int> newMinVertIndex(processors_);
  dynamic_array<int> totVperProc(processors_);
  dynamic_array<int> minNewIndexOnProc(processors_);

  dynamic_array<int> idxIntoSendArray(processors_);
  dynamic_array<int> copyOfReq;

  for (i = 0; i < processors_; ++i)
    send_lens_[i] = 0;

  /*
    compute the prefix sums for vertices going to different processors
    use prefix sum to determine the new vertex indices for local vertices
    newIndex[v] = newVertIndex[vToProc[v]] + minIndexOnProc[vToProc[v]]
  */

  MPI_Allreduce(localVPerProc, totVperProc.data(), processors_, MPI_INT,
                MPI_SUM, comm);
  MPI_Scan(localVPerProc, newMinVertIndex.data(), processors_, MPI_INT,
           MPI_SUM, comm);

  for (i = 0; i < processors_; ++i)
    newMinVertIndex[i] -= localVPerProc[i];

#ifdef DEBUG_HYPERGRAPH
  if (rank_ == 0)
    for (i = 0; i < processors_; ++i)
      assert(newMinVertIndex[i] == 0);
#endif

  dynamic_array<int> maxNewIndexOnProc(processors_);

  minNewIndexOnProc[0] = 0;
  maxNewIndexOnProc[0] = totVperProc[0];

  for (i = 1; i < processors_; ++i) {
    minNewIndexOnProc[i] = minNewIndexOnProc[i - 1] + totVperProc[i - 1];
    maxNewIndexOnProc[i] = minNewIndexOnProc[i] + totVperProc[i];
  }

  for (i = 0; i < processors_; ++i)
    newMinVertIndex[i] = newMinVertIndex[i] + minNewIndexOnProc[i];

  for (i = 0; i < numLocalVertices; ++i) {
    oldIndexToNew[i] = newMinVertIndex[vToProc[i]]++;
#ifdef DEBUG_HYPERGRAPH
    assert(oldIndexToNew[i] >= 0 && oldIndexToNew[i] < numTotalVertices);
#endif
  }

  ds::bit_field sentRequests(numTotalVertices);
  sentRequests.unset();

  /* compute vertices required to transform pinlist */

  for (i = 0; i < numLocalPins; ++i) {
    vertex = localPins[i];

    if (vertex < minVertexIndex || vertex >= maxLocalVertex) {
      if (!sentRequests(vertex)) {
        j = std::min(vertex / vFinePerProc, processors_ - 1);
        data_out_sets_[j]->assign(send_lens_[j]++, vertex);
        sentRequests.set(vertex);
      }
    }
  }

  /* compute number of elements to send to other procs */

  j = 0;
  for (i = 0; i < processors_; ++i) {
    send_displs_[i] = j;
    j += send_lens_[i];
  }

  send_array_.reserve(j);
  copyOfReq.reserve(j);
  totalToSend = j;

  j = 0;
  for (ij = 0; ij < processors_; ++ij) {
    array = data_out_sets_[ij]->data();
    arrayLen = send_lens_[ij];

    for (i = 0; i < arrayLen; ++i) {
      vertex = array[i];
      send_array_[j] = vertex;
      copyOfReq[j++] = vertex;
    }
  }

#ifdef DEBUG_HYPERGRAPH
  assert(j == totalToSend);
#endif

  MPI_Alltoall(send_lens_.data(), 1, MPI_INT, receive_lens_.data(), 1, MPI_INT,
               comm);

  /* compute number of elements to receive from other procs */

  j = 0;
  for (i = 0; i < processors_; ++i) {
    receive_displs_[i] = j;
    j += receive_lens_[i];
  }

  receive_array_.reserve(j);
  totalToRecv = j;

  MPI_Alltoallv(send_array_.data(), send_lens_.data(),
                send_displs_.data(), MPI_INT, receive_array_.data(),
                receive_lens_.data(), receive_displs_.data(), MPI_INT, comm);

  /*
    now have received all requests and sent out our requests
    the reply communication will have the dual dimensions
  */

  send_array_.reserve(totalToRecv);

  for (i = 0; i < totalToRecv; ++i) {
#ifdef DEBUG_HYPERGRAPH
    assert(receive_array_[i] >= minVertexIndex &&
           receive_array_[i] < maxLocalVertex);
#endif
    send_array_[i] = oldIndexToNew[receive_array_[i] - minVertexIndex];
#ifdef DEBUG_HYPERGRAPH
    assert(send_array_[i] >= 0 && send_array_[i] < numTotalVertices);
#endif
  }

  receive_array_.reserve(totalToSend);

  MPI_Alltoallv(send_array_.data(), receive_lens_.data(),
                receive_displs_.data(), MPI_INT, receive_array_.data(),
                send_lens_.data(), send_displs_.data(), MPI_INT, comm);

  /*
   now the requested vertices are in the copyOfReq data_
   while their corresponding matchVector values are in
   the corresponding location in the receive_array_
  */

  if (numLocalPins < numTotalVertices / 2) {
    ds::map_from_pos_int<int> storedRequests(numLocalPins);

    for (i = 0; i < totalToSend; ++i) {
#ifdef DEBUG_HYPERGRAPH
      assert(receive_array_[i] >= 0 && receive_array_[i] < numTotalVertices);
#endif
      storedRequests.insert(copyOfReq[i], receive_array_[i]);
    }

#ifdef DEBUG_HYPERGRAPH
    for (i = 0; i < numLocalVertices; ++i)
      assert(oldIndexToNew[i] >= 0 && oldIndexToNew[i] < numTotalVertices);
#endif

    for (i = 0; i < numLocalPins; ++i) {
      vertex = localPins[i];

      if (vertex >= minVertexIndex && vertex < maxLocalVertex) {
        localPins[i] = oldIndexToNew[vertex - minVertexIndex];
      } else {
        localPins[i] = storedRequests.get(vertex);
#ifdef DEBUG_HYPERGRAPH
        assert(localPins[i] >= 0 && localPins[i] < numTotalVertices);
#endif
      }
    }
  } else {
    dynamic_array<int> nonLocalMatches(numTotalVertices);

    for (i = 0; i < totalToSend; ++i) {
#ifdef DEBUG_HYPERGRAPH
      assert(receive_array_[i] >= 0 && receive_array_[i] < numTotalVertices);
#endif
      nonLocalMatches[copyOfReq[i]] = receive_array_[i];
    }

#ifdef DEBUG_HYPERGRAPH
    for (i = 0; i < numLocalVertices; ++i)
      assert(oldIndexToNew[i] >= 0 && oldIndexToNew[i] < numTotalVertices);
#endif

    for (i = 0; i < numLocalPins; ++i) {
      vertex = localPins[i];

      if (vertex >= minVertexIndex && vertex < maxLocalVertex) {
        localPins[i] = oldIndexToNew[vertex - minVertexIndex];
      } else {
        localPins[i] = nonLocalMatches[vertex];
#ifdef DEBUG_HYPERGRAPH
        assert(localPins[i] >= 0 && localPins[i] < numTotalVertices);
#endif
      }
    }
  }

  /*
    now actually move the vertices to
    their new destination processor
  */

  /* compute the send displacements */

  j = 0;
  i = 0;

  if (numPartitions == 0) {
    totalToSend = Shiftl(numLocalVertices, 1);

    for (; i < processors_; ++i) {
      send_displs_[i] = j;
      idxIntoSendArray[i] = j;
      send_lens_[i] = Shiftl(localVPerProc[i], 1);
      j += send_lens_[i];
    }
  } else {
    totalToSend = numLocalVertices * 3;

    for (; i < processors_; ++i) {
      send_displs_[i] = j;
      idxIntoSendArray[i] = j;
      send_lens_[i] = localVPerProc[i] * 3;
      j += send_lens_[i];
    }
  }

  send_array_.reserve(totalToSend);

  /* compute the send data_ */

  if (numPartitions == 0) {
    for (i = 0; i < numLocalVertices; ++i) {
      j = vToProc[i];
      startOffset = idxIntoSendArray[j];
      send_array_[startOffset] = minVertexIndex + i;
      send_array_[startOffset + 1] = vWeight[i];
      idxIntoSendArray[j] += 2;
    }
  } else {
    for (i = 0; i < numLocalVertices; ++i) {
      j = vToProc[i];
      startOffset = idxIntoSendArray[j];
      send_array_[startOffset] = minVertexIndex + i;
      send_array_[startOffset + 1] = vWeight[i];
      send_array_[startOffset + 2] = partitionVector[i];
      idxIntoSendArray[j] += 3;
    }
  }

  /*
    compute communication dimensions
    and carry out communication
  */

  MPI_Alltoall(send_lens_.data(), 1, MPI_INT, receive_lens_.data(), 1, MPI_INT,
               comm);

  j = 0;
  for (i = 0; i < processors_; ++i) {
    receive_displs_[i] = j;
    j += receive_lens_[i];
  }

#ifdef DEBUG_HYPERGRAPH
  if (numPartitions == 0)
    assert(j == Shiftl(totVperProc[rank_], 1));
  else
    assert(j == totVperProc[rank_] * 3);
#endif

  receive_array_.reserve(j);
  totalToRecv = j;

  MPI_Alltoallv(send_array_.data(), send_lens_.data(),
                send_displs_.data(), MPI_INT, receive_array_.data(),
                receive_lens_.data(), receive_displs_.data(), MPI_INT, comm);

  /* change the local vertex information */

  numLocalVertices = totVperProc[rank_];
  minVertexIndex = 0;

  for (i = 0; i < rank_; ++i)
    minVertexIndex += totVperProc[i];

  vToOrigV.reserve(numLocalVertices);
  vWeight.reserve(numLocalVertices);
  matchVector.reserve(numLocalVertices);

  for (i = 0; i < numLocalVertices; ++i)
    matchVector[i] = -1;

  if (numPartitions > 0) {
    partitionVector.reserve(numLocalVertices);
    partitionOffsetsVector[1] = numLocalVertices;
  }

  i = 0;
  j = 0;
  localVertexWt = 0;

  if (numPartitions == 0) {
    while (i < totalToRecv) {
      vToOrigV[j] = receive_array_[i++];
      vWeight[j++] = receive_array_[i];
      localVertexWt += receive_array_[i++];
    }
  } else {
    while (i < totalToRecv) {
      vToOrigV[j] = receive_array_[i++];
      vWeight[j] = receive_array_[i++];
      partitionVector[j] = receive_array_[i++];
      localVertexWt += vWeight[j++];
    }
  }
}

void ParaHypergraph::shuffleVerticesAftRandom(int *vToProc, int *localVPerProc,
                                              int *mapToOrigV, MPI_Comm comm) {
  /* assume that same number of vertices will remain on each processor */

  int i;
  int j;
  int ij;

  int maxLocalVertex = minVertexIndex + numLocalVertices;
  int vFinePerProc = numTotalVertices / processors_;
  int vToOrigVexist = vToOrigV.capacity();
  int totalToSend;
  int totalToRecv;
  int vertex;
  int startOffset;
  int arrayLen;
  int *array;

  dynamic_array<int> oldIndexToNew(numLocalVertices);
  dynamic_array<int> minIndexOfMyVerts(processors_);
  dynamic_array<int> minIndexOnProc(processors_);
  dynamic_array<int> idxIntoSendArray(processors_);
  dynamic_array<int> copyOfReq;

  for (i = 0; i < processors_; ++i)
    send_lens_[i] = 0;

  /*
    compute the prefix sums for vertices going to different processors
    use prefix sum to determine the new vertex indices for local vertices
    newIndex[v] = newVertIndex[vToProc[v]] + minIndexOnProc[vToProc[v]]
  */

  MPI_Scan(localVPerProc, minIndexOfMyVerts.data(), processors_, MPI_INT,
           MPI_SUM, comm);

  for (i = 0; i < processors_; ++i)
    minIndexOfMyVerts[i] -= localVPerProc[i];

  j = 0;
  ij = numTotalVertices / processors_;

  for (i = 0; i < processors_; ++i) {
    minIndexOnProc[i] = j;
    j += ij;
  }

  for (i = 0; i < processors_; ++i)
    minIndexOfMyVerts[i] += minIndexOnProc[i];

  for (i = 0; i < numLocalVertices; ++i) {
    oldIndexToNew[i] = minIndexOfMyVerts[vToProc[i]]++;
#ifdef DEBUG_HYPERGRAPH
    assert(oldIndexToNew[i] >= 0 && oldIndexToNew[i] < numTotalVertices);
#endif
  }

  ds::bit_field sentRequests(numTotalVertices);
  sentRequests.unset();

  /* compute all the requests for new indices of remote vertices */

  for (i = 0; i < numLocalPins; ++i) {
    vertex = localPins[i];

    if (vertex < minVertexIndex || vertex >= maxLocalVertex) {
      if (!sentRequests(vertex)) {
        j = std::min(vertex / vFinePerProc, processors_ - 1);
        data_out_sets_[j]->assign(send_lens_[j]++, vertex);
        sentRequests.set(vertex);
      }
    }
  }

  /* compute number of elements to send to other procs */

  j = 0;
  for (i = 0; i < processors_; ++i) {
    send_displs_[i] = j;
    j += send_lens_[i];
  }

  send_array_.reserve(j);
  copyOfReq.reserve(j);
  totalToSend = j;

  j = 0;
  for (ij = 0; ij < processors_; ++ij) {
    array = data_out_sets_[ij]->data();
    arrayLen = send_lens_[ij];

    for (i = 0; i < arrayLen; ++i) {
      vertex = array[i];
      send_array_[j] = vertex;
      copyOfReq[j++] = vertex;
    }
  }

#ifdef DEBUG_HYPERGRAPH
  assert(j == totalToSend);
#endif

  MPI_Alltoall(send_lens_.data(), 1, MPI_INT, receive_lens_.data(), 1, MPI_INT,
               comm);

  /* compute number of elements to receive from other procs */

  j = 0;
  for (i = 0; i < processors_; ++i) {
    receive_displs_[i] = j;
    j += receive_lens_[i];
  }

  receive_array_.reserve(j);
  totalToRecv = j;

  MPI_Alltoallv(send_array_.data(), send_lens_.data(),
                send_displs_.data(), MPI_INT, receive_array_.data(),
                receive_lens_.data(), receive_displs_.data(), MPI_INT, comm);

  /*
    now have received all requests and sent out our requests
    the reply communication will have the dual dimensions
  */

  send_array_.reserve(totalToRecv);

  for (i = 0; i < totalToRecv; ++i) {
#ifdef DEBUG_HYPERGRAPH
    assert(receive_array_[i] >= minVertexIndex &&
           receive_array_[i] < maxLocalVertex);
#endif
    send_array_[i] = oldIndexToNew[receive_array_[i] - minVertexIndex];
#ifdef DEBUG_HYPERGRAPH
    assert(send_array_[i] >= 0 && send_array_[i] < numTotalVertices);
#endif
  }

  receive_array_.reserve(totalToSend);

  MPI_Alltoallv(send_array_.data(), receive_lens_.data(),
                receive_displs_.data(), MPI_INT, receive_array_.data(),
                send_lens_.data(), send_displs_.data(), MPI_INT, comm);

  /*
     now the requested vertices are in the copyOfReq data_
     while their corresponding matchVector values are in
     the corresponding location in the receive_array_
  */

  if (numLocalPins < numTotalVertices / 2) {
    ds::map_from_pos_int<int> storedRequests(numLocalPins);

    for (i = 0; i < totalToSend; ++i) {
#ifdef DEBUG_HYPERGRAPH
      assert(receive_array_[i] >= 0 && receive_array_[i] < numTotalVertices);
#endif
      storedRequests.insert(copyOfReq[i], receive_array_[i]);
    }

#ifdef DEBUG_HYPERGRAPH
    for (i = 0; i < numLocalVertices; ++i)
      assert(oldIndexToNew[i] >= 0 && oldIndexToNew[i] < numTotalVertices);
#endif

    /*  now we will convert the local pin list and hash keys */

    for (i = 0; i < numLocalPins; ++i) {
      vertex = localPins[i];

      if (vertex >= minVertexIndex && vertex < maxLocalVertex) {
        localPins[i] = oldIndexToNew[vertex - minVertexIndex];
      } else {
        localPins[i] = storedRequests.get(vertex);
#ifdef DEBUG_HYPERGRAPH
        assert(localPins[i] >= 0 && localPins[i] < numTotalVertices);
#endif
      }
    }
  } else {
    dynamic_array<int> storedRequests(numTotalVertices);

    for (i = 0; i < totalToSend; ++i) {
#ifdef DEBUG_HYPERGRAPH
      assert(receive_array_[i] >= 0 && receive_array_[i] < numTotalVertices);
#endif
      storedRequests[copyOfReq[i]] = receive_array_[i];
    }

#ifdef DEBUG_HYPERGRAPH
    for (i = 0; i < numLocalVertices; ++i)
      assert(oldIndexToNew[i] >= 0 && oldIndexToNew[i] < numTotalVertices);
#endif

    /*  now we will convert the local pin list and hash keys */

    for (i = 0; i < numLocalPins; ++i) {
      vertex = localPins[i];

      if (vertex >= minVertexIndex && vertex < maxLocalVertex) {
        localPins[i] = oldIndexToNew[vertex - minVertexIndex];
      } else {
        localPins[i] = storedRequests[vertex];
#ifdef DEBUG_HYPERGRAPH
        assert(localPins[i] >= 0 && localPins[i] < numTotalVertices);
#endif
      }
    }
  }

  /*
    now actually move the vertices to
    their new destination processor
  */

  /* compute the send displacements */

  j = 0;

  if (!vToOrigVexist) {
    ij = numPartitions + 2;
    totalToSend = numLocalVertices * ij;
#ifdef DEBUG_HYPERGRAPH
    assert(totalToSend >= 0);
#endif
    for (i = 0; i < processors_; ++i) {
      send_displs_[i] = j;
      idxIntoSendArray[i] = j;
      send_lens_[i] = localVPerProc[i] * ij;
      j += send_lens_[i];
    }
  } else {
    ij = numPartitions + 3;
    totalToSend = numLocalVertices * ij;
#ifdef DEBUG_HYPERGRAPH
    assert(totalToSend >= 0);
#endif
    for (i = 0; i < processors_; ++i) {
      send_displs_[i] = j;
      idxIntoSendArray[i] = j;
      send_lens_[i] = localVPerProc[i] * ij;
      j += send_lens_[i];
    }
  }

#ifdef DEBUG_HYPERGRAPH
  assert(j == totalToSend);
  assert(totalToSend >= 0);
#endif

  send_array_.reserve(totalToSend);

  /* compute the send data_ */

  if (!vToOrigVexist) {
    for (i = 0; i < numLocalVertices; ++i) {
      j = vToProc[i];
      startOffset = idxIntoSendArray[j];
      send_array_[startOffset++] = vWeight[i];
      send_array_[startOffset++] = mapToOrigV[i];

      for (ij = 0; ij < numPartitions; ++ij)
        send_array_[startOffset++] =
            partitionVector[partitionOffsetsVector[ij] + i];

      idxIntoSendArray[j] += (ij + 2);
    }
  } else {
    for (i = 0; i < numLocalVertices; ++i) {
      j = vToProc[i];
      startOffset = idxIntoSendArray[j];
      send_array_[startOffset++] = vToOrigV[i];
      send_array_[startOffset++] = vWeight[i];
      send_array_[startOffset++] = mapToOrigV[i];

      for (ij = 0; ij < numPartitions; ++ij)
        send_array_[startOffset++] =
            partitionVector[partitionOffsetsVector[ij] + i];

      idxIntoSendArray[j] += (ij + 3);
    }
  }

  /*
    compute communication dimensions
    and carry out communication
  */

  MPI_Alltoall(send_lens_.data(), 1, MPI_INT, receive_lens_.data(), 1, MPI_INT,
               comm);

  j = 0;
  for (i = 0; i < processors_; ++i) {
    receive_displs_[i] = j;
    j += receive_lens_[i];
  }

#ifdef DEBUG_HYPERGRAPH
  if (!vToOrigVexist)
    assert(j == numLocalVertices * (numPartitions + 2));
  else
    assert(j == numLocalVertices * (numPartitions + 3));
#endif

  receive_array_.reserve(j);
  totalToRecv = j;

  MPI_Alltoallv(send_array_.data(), send_lens_.data(),
                send_displs_.data(), MPI_INT, receive_array_.data(),
                receive_lens_.data(), receive_displs_.data(), MPI_INT, comm);

  /* change the local vertex information */

  i = 0;
  j = 0;
  localVertexWt = 0;

  if (!vToOrigVexist) {
    while (i < totalToRecv) {
      vWeight[j] = receive_array_[i];
      localVertexWt += receive_array_[i++];
      mapToOrigV[j] = receive_array_[i++];

      for (ij = 0; ij < numPartitions; ++ij)
        partitionVector[partitionOffsetsVector[ij] + j] = receive_array_[i++];

      ++j;
    }
  } else {
    vToOrigV.reserve(numLocalVertices);

    while (i < totalToRecv) {
      vToOrigV[j] = receive_array_[i++];
      vWeight[j] = receive_array_[i];
      localVertexWt += receive_array_[i++];
      mapToOrigV[j] = receive_array_[i++];

      for (ij = 0; ij < numPartitions; ++ij)
        partitionVector[partitionOffsetsVector[ij] + j] = receive_array_[i++];

      ++j;
    }
  }
}

void ParaHypergraph::shuffleVerticesAftRandom(int *vToProc, int *localVPerProc,
                                              int *mapToInterV, int *mapToOrigV,
                                              MPI_Comm comm) {
  /* assume that same number of vertices will remain on each processor */

  int i;
  int j;
  int ij;

  int maxLocalVertex = minVertexIndex + numLocalVertices;
  int vFinePerProc = numTotalVertices / processors_;
  int vToOrigVexist = vToOrigV.capacity();
  int totalToSend;
  int totalToRecv;
  int vertex;
  int startOffset;
  int arrayLen;
  int *array;

  dynamic_array<int> oldIndexToNew(numLocalVertices);
  dynamic_array<int> minIndexOfMyVerts(processors_);
  dynamic_array<int> minIndexOnProc(processors_);
  dynamic_array<int> idxIntoSendArray(processors_);
  dynamic_array<int> copyOfReq;

  for (i = 0; i < processors_; ++i)
    send_lens_[i] = 0;

  /*
    compute the prefix sums for vertices going to different processors
    use prefix sum to determine the new vertex indices for local vertices
    newIndex[v] = newVertIndex[vToProc[v]] + minIndexOnProc[vToProc[v]]
  */

  MPI_Scan(localVPerProc, minIndexOfMyVerts.data(), processors_, MPI_INT,
           MPI_SUM, comm);

  for (i = 0; i < processors_; ++i)
    minIndexOfMyVerts[i] -= localVPerProc[i];

  j = 0;
  ij = numTotalVertices / processors_;

  for (i = 0; i < processors_; ++i) {
    minIndexOnProc[i] = j;
    j += ij;
  }

  for (i = 0; i < processors_; ++i)
    minIndexOfMyVerts[i] += minIndexOnProc[i];

  for (i = 0; i < numLocalVertices; ++i) {
    oldIndexToNew[i] = minIndexOfMyVerts[vToProc[i]]++;
#ifdef DEBUG_HYPERGRAPH
    assert(oldIndexToNew[i] >= 0 && oldIndexToNew[i] < numTotalVertices);
#endif
  }

  /* now need to convert the pinlist and the hyperedge hash keys */

  ds::bit_field sentRequests(numTotalVertices);
  sentRequests.unset();

  /* compute all the requests for new indices of remote vertices  */

  for (i = 0; i < numLocalPins; ++i) {
    vertex = localPins[i];

    if (vertex < minVertexIndex || vertex >= maxLocalVertex) {
      if (!sentRequests(vertex)) {
        j = std::min(vertex / vFinePerProc, processors_ - 1);
        data_out_sets_[j]->assign(send_lens_[j]++, vertex);
        sentRequests.set(vertex);
      }
    }
  }

  /* compute number of elements to send to other procs */

  j = 0;
  for (i = 0; i < processors_; ++i) {
    send_displs_[i] = j;
    j += send_lens_[i];
  }

  send_array_.reserve(j);
  copyOfReq.reserve(j);
  totalToSend = j;

  j = 0;
  for (ij = 0; ij < processors_; ++ij) {
    array = data_out_sets_[ij]->data();
    arrayLen = send_lens_[ij];

    for (i = 0; i < arrayLen; ++i) {
      vertex = array[i];
      send_array_[j] = vertex;
      copyOfReq[j++] = vertex;
    }
  }

#ifdef DEBUG_HYPERGRAPH
  assert(j == totalToSend);
#endif

  MPI_Alltoall(send_lens_.data(), 1, MPI_INT, receive_lens_.data(), 1, MPI_INT,
               comm);

  /* compute number of elements to receive from other procs */

  j = 0;
  for (i = 0; i < processors_; ++i) {
    receive_displs_[i] = j;
    j += receive_lens_[i];
  }

  receive_array_.reserve(j);
  totalToRecv = j;

  MPI_Alltoallv(send_array_.data(), send_lens_.data(),
                send_displs_.data(), MPI_INT, receive_array_.data(),
                receive_lens_.data(), receive_displs_.data(), MPI_INT, comm);

  /*
    now have received all requests and sent out our requests
    the reply communication will have the dual dimensions
  */

  send_array_.reserve(totalToRecv);

  for (i = 0; i < totalToRecv; ++i) {
#ifdef DEBUG_HYPERGRAPH
    assert(receive_array_[i] >= minVertexIndex &&
           receive_array_[i] < maxLocalVertex);
#endif
    send_array_[i] = oldIndexToNew[receive_array_[i] - minVertexIndex];
#ifdef DEBUG_HYPERGRAPH
    assert(send_array_[i] >= 0 && send_array_[i] < numTotalVertices);
#endif
  }

  receive_array_.reserve(totalToSend);

  MPI_Alltoallv(send_array_.data(), receive_lens_.data(),
                receive_displs_.data(), MPI_INT, receive_array_.data(),
                send_lens_.data(), send_displs_.data(), MPI_INT, comm);

  /*
  now the requested vertices are in the copyOfReq data_
   while their corresponding matchVector values are in
   the corresponding location in the receive_array_
  */

  if (numLocalPins < numTotalVertices / 2) {
    ds::map_from_pos_int<int> storedRequests(numLocalPins);

    for (i = 0; i < totalToSend; ++i) {
#ifdef DEBUG_HYPERGRAPH
      assert(receive_array_[i] >= 0 && receive_array_[i] < numTotalVertices);
#endif
      storedRequests.insert(copyOfReq[i], receive_array_[i]);
    }

/* now we will convert the local pin list and hash keys */

#ifdef DEBUG_HYPERGRAPH
    for (i = 0; i < numLocalVertices; ++i)
      assert(oldIndexToNew[i] >= 0 && oldIndexToNew[i] < numTotalVertices);
#endif

    for (i = 0; i < numLocalPins; ++i) {
      vertex = localPins[i];

      if (vertex >= minVertexIndex && vertex < maxLocalVertex) {
        localPins[i] = oldIndexToNew[vertex - minVertexIndex];
      } else {
        localPins[i] = storedRequests.get(vertex);
#ifdef DEBUG_HYPERGRAPH
        assert(localPins[i] >= 0 && localPins[i] < numTotalVertices);
#endif
      }
    }
  } else {
    dynamic_array<int> storedRequests(numTotalVertices);

    for (i = 0; i < totalToSend; ++i) {
#ifdef DEBUG_HYPERGRAPH
      assert(receive_array_[i] >= 0 && receive_array_[i] < numTotalVertices);
#endif
      storedRequests[copyOfReq[i]] = receive_array_[i];
    }

/* now we will convert the local pin list and hash keys */

#ifdef DEBUG_HYPERGRAPH
    for (i = 0; i < numLocalVertices; ++i)
      assert(oldIndexToNew[i] >= 0 && oldIndexToNew[i] < numTotalVertices);
#endif

    for (i = 0; i < numLocalPins; ++i) {
      vertex = localPins[i];

      if (vertex >= minVertexIndex && vertex < maxLocalVertex) {
        localPins[i] = oldIndexToNew[vertex - minVertexIndex];
      } else {
        localPins[i] = storedRequests[vertex];
#ifdef DEBUG_HYPERGRAPH
        assert(localPins[i] >= 0 && localPins[i] < numTotalVertices);
#endif
      }
    }
  }

  /*
    now actually move the vertices to
    their new destination processor
  */

  /* compute the send displacements */

  j = 0;

  if (!vToOrigVexist) {
    ij = numPartitions + 3;
    totalToSend = numLocalVertices * ij;

    for (i = 0; i < processors_; ++i) {
      send_displs_[i] = j;
      idxIntoSendArray[i] = j;
      send_lens_[i] = localVPerProc[i] * ij;
      j += send_lens_[i];
    }
  } else {
    ij = numPartitions + 4;
    totalToSend = numLocalVertices * ij;

    for (i = 0; i < processors_; ++i) {
      send_displs_[i] = j;
      idxIntoSendArray[i] = j;
      send_lens_[i] = localVPerProc[i] * ij;
      j += send_lens_[i];
    }
  }

  send_array_.reserve(totalToSend);

  /* compute the send data_ */

  if (!vToOrigVexist) {
    for (i = 0; i < numLocalVertices; ++i) {
      j = vToProc[i];
      startOffset = idxIntoSendArray[j];
      send_array_[startOffset++] = vWeight[i];
      send_array_[startOffset++] = mapToInterV[i];
      send_array_[startOffset++] = mapToOrigV[i];

      for (ij = 0; ij < numPartitions; ++ij)
        send_array_[startOffset++] =
            partitionVector[partitionOffsetsVector[ij] + i];

      idxIntoSendArray[j] += (3 + ij);
    }
  } else {
    for (i = 0; i < numLocalVertices; ++i) {
      j = vToProc[i];
      startOffset = idxIntoSendArray[j];
      send_array_[startOffset++] = vToOrigV[i];
      send_array_[startOffset++] = vWeight[i];
      send_array_[startOffset++] = mapToInterV[i];
      send_array_[startOffset++] = mapToOrigV[i];

      for (ij = 0; ij < numPartitions; ++ij)
        send_array_[startOffset++] =
            partitionVector[partitionOffsetsVector[ij] + i];

      idxIntoSendArray[j] += (4 + ij);
    }
  }

  /*
     compute communication dimensions
     and carry out communication
  */

  MPI_Alltoall(send_lens_.data(), 1, MPI_INT, receive_lens_.data(), 1, MPI_INT,
               comm);

  j = 0;
  for (i = 0; i < processors_; ++i) {
    receive_displs_[i] = j;
    j += receive_lens_[i];
  }

#ifdef DEBUG_HYPERGRAPH
  if (!vToOrigVexist)
    assert(j == numLocalVertices * (3 + numPartitions));
  if (vToOrigVexist)
    assert(j == numLocalVertices * (4 + numPartitions));
#endif

  receive_array_.reserve(j);
  totalToRecv = j;

  MPI_Alltoallv(send_array_.data(), send_lens_.data(),
                send_displs_.data(), MPI_INT, receive_array_.data(),
                receive_lens_.data(), receive_displs_.data(), MPI_INT, comm);

  /* change the local vertex information */

  i = 0;
  j = 0;
  localVertexWt = 0;

  if (!vToOrigVexist) {
    while (i < totalToRecv) {
      vWeight[j] = receive_array_[i];
      localVertexWt += receive_array_[i++];
      mapToInterV[j] = receive_array_[i++];
      mapToOrigV[j] = receive_array_[i++];

      for (ij = 0; ij < numPartitions; ++ij)
        partitionVector[partitionOffsetsVector[ij] + j] = receive_array_[i++];

      ++j;
    }
  } else {
    vToOrigV.reserve(numLocalVertices);

    while (i < totalToRecv) {
      vToOrigV[j] = receive_array_[i++];
      vWeight[j] = receive_array_[i];
      localVertexWt += receive_array_[i++];
      mapToInterV[j] = receive_array_[i++];
      mapToOrigV[j] = receive_array_[i++];

      for (ij = 0; ij < numPartitions; ++ij)
        partitionVector[partitionOffsetsVector[ij] + j] = receive_array_[i++];

      ++j;
    }
  }
}

void ParaHypergraph::shuffleVerticesAftRandom(int *vToProc, int *localVPerProc,
                                              ParaHypergraph &fineG,
                                              MPI_Comm comm) {
#ifdef DEBUG_HYPERGRAPH
  assert(numPartitions > 0);
#endif

  int i;
  int j;
  int ij;

  int maxLocalVertex = minVertexIndex + numLocalVertices;
  int vFinePerProc = numTotalVertices / processors_;
  int totalToSend;
  int totalToRecv;
  int vertex;
  int startOffset;
  int arrayLen;
  int vToOrigVexist = vToOrigV.capacity();
  int *array;

  dynamic_array<int> oldIndexToNew(numLocalVertices);
  dynamic_array<int> minIndexOfMyVerts(processors_);
  dynamic_array<int> minIndexOnProc(processors_);
  dynamic_array<int> idxIntoSendArray(processors_);
  dynamic_array<int> copyOfReq;

  /* finer graph structs */

  int numLocalFineVertices = fineG.getNumLocalVertices();
  int *fineMatchVector = fineG.getMatchVectorArray();

  for (i = 0; i < processors_; ++i)
    send_lens_[i] = 0;

  MPI_Scan(localVPerProc, minIndexOfMyVerts.data(), processors_, MPI_INT,
           MPI_SUM, comm);

  for (i = 0; i < processors_; ++i)
    minIndexOfMyVerts[i] -= localVPerProc[i];

  j = 0;
  ij = numTotalVertices / processors_;

  for (i = 0; i < processors_; ++i) {
    minIndexOnProc[i] = j;
    j += ij;
  }

  for (i = 0; i < processors_; ++i)
    minIndexOfMyVerts[i] += minIndexOnProc[i];

  for (i = 0; i < numLocalVertices; ++i) {
    oldIndexToNew[i] = minIndexOfMyVerts[vToProc[i]]++;
#ifdef DEBUG_HYPERGRAPH
    assert(oldIndexToNew[i] >= 0 && oldIndexToNew[i] < numTotalVertices);
#endif
  }

  /* now need to convert the pinlist and the hyperedge hash keys */

  ds::bit_field sentRequests(numTotalVertices);
  sentRequests.unset();

  /*
     compute all the requests for new indices of remote vertices
     because we need to modify the match vector of the finer graph,
     this needs to include those vertices that are in this match
     vector as well
  */

  for (i = 0; i < numLocalPins; ++i) {
    vertex = localPins[i];

    if (vertex < minVertexIndex || vertex >= maxLocalVertex) {
      if (!sentRequests(vertex)) {
        j = std::min(vertex / vFinePerProc, processors_ - 1);
        data_out_sets_[j]->assign(send_lens_[j]++, vertex);
        sentRequests.set(vertex);
      }
    }
  }

  for (i = 0; i < numLocalFineVertices; ++i) {
    vertex = fineMatchVector[i];

    if (vertex < minVertexIndex || vertex >= maxLocalVertex) {
      if (!sentRequests(vertex)) {
        j = std::min(vertex / vFinePerProc, processors_ - 1);
        data_out_sets_[j]->assign(send_lens_[j]++, vertex);
        sentRequests.set(vertex);
      }
    }
  }

  /* compute number of elements to send to other procs */

  j = 0;
  for (i = 0; i < processors_; ++i) {
    send_displs_[i] = j;
    j += send_lens_[i];
  }

  send_array_.reserve(j);
  copyOfReq.reserve(j);
  totalToSend = j;

  j = 0;
  for (ij = 0; ij < processors_; ++ij) {
    array = data_out_sets_[ij]->data();
    arrayLen = send_lens_[ij];

    for (i = 0; i < arrayLen; ++i) {
      vertex = array[i];
      send_array_[j] = vertex;
      copyOfReq[j++] = vertex;
    }
  }

#ifdef DEBUG_HYPERGRAPH
  assert(j == totalToSend);
#endif

  MPI_Alltoall(send_lens_.data(), 1, MPI_INT, receive_lens_.data(), 1, MPI_INT,
               comm);

  /* compute number of elements to receive from other procs */

  j = 0;
  for (i = 0; i < processors_; ++i) {
    receive_displs_[i] = j;
    j += receive_lens_[i];
  }

  receive_array_.reserve(j);
  totalToRecv = j;

  MPI_Alltoallv(send_array_.data(), send_lens_.data(),
                send_displs_.data(), MPI_INT, receive_array_.data(),
                receive_lens_.data(), receive_displs_.data(), MPI_INT, comm);

  /*
    now have received all requests and sent out our requests
    the reply communication will have the dual dimensions
  */

  send_array_.reserve(totalToRecv);

  for (i = 0; i < totalToRecv; ++i) {
#ifdef DEBUG_HYPERGRAPH
    assert(receive_array_[i] >= minVertexIndex &&
           receive_array_[i] < maxLocalVertex);
#endif
    send_array_[i] = oldIndexToNew[receive_array_[i] - minVertexIndex];
#ifdef DEBUG_HYPERGRAPH
    assert(send_array_[i] >= 0 && send_array_[i] < numTotalVertices);
#endif
  }

  receive_array_.reserve(totalToSend);

  MPI_Alltoallv(send_array_.data(), receive_lens_.data(),
                receive_displs_.data(), MPI_INT, receive_array_.data(),
                send_lens_.data(), send_displs_.data(), MPI_INT, comm);

  /*
    now the requested vertices are in the copyOfReq data_
    while their corresponding matchVector values are in
  */

  if (numLocalPins < numTotalVertices / 2) {
    ds::map_from_pos_int<int> storedRequests(numLocalPins);

    for (i = 0; i < totalToSend; ++i) {
#ifdef DEBUG_HYPERGRAPH
      assert(receive_array_[i] >= 0 && receive_array_[i] < numTotalVertices);
#endif
      storedRequests.insert(copyOfReq[i], receive_array_[i]);
    }

/* now we will convert the local pin list and hash keys */

#ifdef DEBUG_HYPERGRAPH
    for (i = 0; i < numLocalVertices; ++i)
      assert(oldIndexToNew[i] >= 0 && oldIndexToNew[i] < numTotalVertices);
#endif

    for (i = 0; i < numLocalPins; ++i) {
      vertex = localPins[i];

      if (vertex >= minVertexIndex && vertex < maxLocalVertex) {
        localPins[i] = oldIndexToNew[vertex - minVertexIndex];
      } else {
        localPins[i] = storedRequests.get(vertex);
#ifdef DEBUG_HYPERGRAPH
        assert(localPins[i] >= 0 && localPins[i] < numTotalVertices);
#endif
      }
    }

    /* now we will convert the fine graph match vector */

    for (i = 0; i < numLocalFineVertices; ++i) {
      vertex = fineMatchVector[i];

      if (vertex >= minVertexIndex && vertex < maxLocalVertex) {
        fineMatchVector[i] = oldIndexToNew[vertex - minVertexIndex];
      } else {
        fineMatchVector[i] = storedRequests.get(vertex);
      }
#ifdef DEBUG_HYPERGRAPH
      assert(fineMatchVector[i] >= 0 && fineMatchVector[i] < numTotalVertices);
#endif
    }
  } else {
    dynamic_array<int> storedRequests(numTotalVertices);

    for (i = 0; i < totalToSend; ++i) {
#ifdef DEBUG_HYPERGRAPH
      assert(receive_array_[i] >= 0 && receive_array_[i] < numTotalVertices);
#endif
      storedRequests[copyOfReq[i]] = receive_array_[i];
    }

/* now we will convert the local pin list and hash keys */

#ifdef DEBUG_HYPERGRAPH
    for (i = 0; i < numLocalVertices; ++i)
      assert(oldIndexToNew[i] >= 0 && oldIndexToNew[i] < numTotalVertices);
#endif

    for (i = 0; i < numLocalPins; ++i) {
      vertex = localPins[i];

      if (vertex >= minVertexIndex && vertex < maxLocalVertex) {
        localPins[i] = oldIndexToNew[vertex - minVertexIndex];
      } else {
        localPins[i] = storedRequests[vertex];
#ifdef DEBUG_HYPERGRAPH
        assert(localPins[i] >= 0 && localPins[i] < numTotalVertices);
#endif
      }
    }

    /* now we will convert the fine graph match vector */

    for (i = 0; i < numLocalFineVertices; ++i) {
      vertex = fineMatchVector[i];

      if (vertex >= minVertexIndex && vertex < maxLocalVertex) {
        fineMatchVector[i] = oldIndexToNew[vertex - minVertexIndex];
      } else {
        fineMatchVector[i] = storedRequests[vertex];
      }
#ifdef DEBUG_HYPERGRAPH
      assert(fineMatchVector[i] >= 0 && fineMatchVector[i] < numTotalVertices);
#endif
    }
  }

  /*
    now actually move the vertices to
    their new destination processor
  */

  /* compute the send displacements */

  j = 0;

  if (vToOrigVexist > 0) {
    ij = numPartitions + 2;
    totalToSend = numLocalVertices * ij;

    for (i = 0; i < processors_; ++i) {
      send_displs_[i] = j;
      idxIntoSendArray[i] = j;
      send_lens_[i] = localVPerProc[i] * ij;
      j += send_lens_[i];
    }
  } else {
    ij = numPartitions + 1;
    totalToSend = numLocalVertices * ij;

    for (i = 0; i < processors_; ++i) {
      send_displs_[i] = j;
      idxIntoSendArray[i] = j;
      send_lens_[i] = localVPerProc[i] * ij;
      j += send_lens_[i];
    }
  }

#ifdef DEBUG_HYPERGRAPH
  assert(j == totalToSend);
#endif

  send_array_.reserve(totalToSend);

  /* compute the send data_ */

  if (vToOrigVexist > 0) {
    for (i = 0; i < numLocalVertices; ++i) {
      j = vToProc[i];

      startOffset = idxIntoSendArray[j];
      send_array_[startOffset++] = vToOrigV[i];
      send_array_[startOffset++] = vWeight[i];

      for (ij = 0; ij < numPartitions; ++ij)
        send_array_[startOffset++] =
            partitionVector[partitionOffsetsVector[ij] + i];

      idxIntoSendArray[j] += (ij + 2);
    }
  } else {
    for (i = 0; i < numLocalVertices; ++i) {
      j = vToProc[i];

      startOffset = idxIntoSendArray[j];
      send_array_[startOffset++] = vWeight[i];

      for (ij = 0; ij < numPartitions; ++ij)
        send_array_[startOffset++] =
            partitionVector[partitionOffsetsVector[ij] + i];

      idxIntoSendArray[j] += (ij + 1);
    }
  }

  /*
     compute communication dimensions
     and carry out communication
  */

  MPI_Alltoall(send_lens_.data(), 1, MPI_INT, receive_lens_.data(), 1, MPI_INT,
               comm);

  j = 0;
  for (i = 0; i < processors_; ++i) {
    receive_displs_[i] = j;
    j += receive_lens_[i];
  }

#ifdef DEBUG_HYPERGRAPH
  if (vToOrigVexist > 0)
    assert(j == numLocalVertices * (numPartitions + 2));
  else
    assert(j == numLocalVertices * (numPartitions + 1));
#endif

  receive_array_.reserve(j);
  totalToRecv = j;

  MPI_Alltoallv(send_array_.data(), send_lens_.data(),
                send_displs_.data(), MPI_INT, receive_array_.data(),
                receive_lens_.data(), receive_displs_.data(), MPI_INT, comm);

  /* change the local vertex information */

  i = 0;
  j = 0;
  localVertexWt = 0;

  if (vToOrigVexist > 0) {
    while (i < totalToRecv) {
      vToOrigV[j] = receive_array_[i++];
      vWeight[j] = receive_array_[i];
      localVertexWt += receive_array_[i++];

      for (ij = 0; ij < numPartitions; ++ij)
        partitionVector[partitionOffsetsVector[ij] + j] = receive_array_[i++];

      ++j;
    }
  } else {
    while (i < totalToRecv) {
      vWeight[j] = receive_array_[i];
      localVertexWt += receive_array_[i++];

      for (ij = 0; ij < numPartitions; ++ij)
        partitionVector[partitionOffsetsVector[ij] + j] = receive_array_[i++];

      ++j;
    }
  }
}

void ParaHypergraph::shiftVerticesToBalance(MPI_Comm comm) {
#ifdef DEBUG_HYPERGRAPH
  assert(numPartitions == 0);
#endif

  int i;
  int j;

  int vPerProc = numTotalVertices / processors_;
  int maxVertexIndex = minVertexIndex + numLocalVertices;
  int numMyVertices;

  if (rank_ != processors_ - 1)
    numMyVertices = vPerProc;
  else
    numMyVertices = vPerProc + Mod(numTotalVertices, processors_);

  dynamic_array<int> minNewIndex(processors_);
  dynamic_array<int> maxNewIndex(processors_);

  for (i = 0; i < processors_; ++i) {
    if (i == 0) {
      minNewIndex[i] = 0;
      maxNewIndex[i] = vPerProc;
    } else {
      minNewIndex[i] = maxNewIndex[i - 1];
      if (i == processors_ - 1)
        maxNewIndex[i] = numTotalVertices;
      else
        maxNewIndex[i] = minNewIndex[i] + vPerProc;
    }
  }

  /*
    here distinguish cases where VtoOrigV data_ needs
    to be maintained
  */

  if (vToOrigV.capacity() > 0) {
    j = 0;
    send_array_.reserve(numLocalVertices * 3);

    for (i = 0; i < numLocalVertices; ++i) {
      send_array_[j++] = vWeight[i];
      send_array_[j++] = matchVector[i];
      send_array_[j++] = vToOrigV[i];
    }
#ifdef DEBUG_HYPERGRAPH
    assert(j == numLocalVertices * 3);
#endif

    for (i = 0; i < processors_; ++i) {
      if (i == 0)
        send_displs_[i] = 0;
      else
        send_displs_[i] = send_displs_[i - 1] + send_lens_[i - 1];

      send_lens_[i] =
          (std::max(numLocalVertices - (std::max(maxVertexIndex - maxNewIndex[i], 0) +
                                   std::max(minNewIndex[i] - minVertexIndex, 0)),
               0)) *
          3;
    }
  } else {
    j = 0;
    send_array_.reserve(Shiftl(numLocalVertices, 1));

    for (i = 0; i < numLocalVertices; ++i) {
      send_array_[j++] = vWeight[i];
      send_array_[j++] = matchVector[i];
    }
#ifdef DEBUG_HYPERGRAPH
    assert(j == Shiftl(numLocalVertices, 1));
#endif

    for (i = 0; i < processors_; ++i) {
      if (i == 0)
        send_displs_[i] = 0;
      else
        send_displs_[i] = send_displs_[i - 1] + send_lens_[i - 1];

      send_lens_[i] = Shiftl(
          std::max(numLocalVertices - (std::max(maxVertexIndex - maxNewIndex[i], 0) +
                                  std::max(minNewIndex[i] - minVertexIndex, 0)),
              0),
          1);
    }
  }

  MPI_Alltoall(send_lens_.data(), 1, MPI_INT, receive_lens_.data(), 1, MPI_INT,
               comm);

  j = 0;
  for (i = 0; i < processors_; ++i) {
    receive_displs_[i] = j;
    j += receive_lens_[i];
  }

#ifdef DEBUG_HYPERGRAPH
  if (vToOrigV.getLength() > 0)
    assert(j == numMyVertices * 3);
  else
    assert(j == Shiftl(numMyVertices, 1));
#endif

  receive_array_.reserve(j);

  MPI_Alltoallv(send_array_.data(), send_lens_.data(),
                send_displs_.data(), MPI_INT, receive_array_.data(),
                receive_lens_.data(), receive_displs_.data(), MPI_INT, comm);

  numLocalVertices = numMyVertices;
  minVertexIndex = minNewIndex[rank_];

  vWeight.reserve(numLocalVertices);
  matchVector.reserve(numLocalVertices);

  if (vToOrigV.capacity() > 0) {
    vToOrigV.reserve(numLocalVertices);

    j = 0;
    for (i = 0; i < numLocalVertices; ++i) {
      vWeight[i] = receive_array_[j++];
      matchVector[i] = receive_array_[j++];
      vToOrigV[i] = receive_array_[j++];
    }
  } else {
    j = 0;
    for (i = 0; i < numLocalVertices; ++i) {
      vWeight[i] = receive_array_[j++];
      matchVector[i] = receive_array_[j++];
    }
  }
}

int ParaHypergraph::calcCutsize(int numParts, int pNum, MPI_Comm comm) {
#ifdef DEBUG_HYPERGRAPH
  assert(numPartitions > 0);
  assert(pNum >= 0 && pNum < numPartitions);
#endif

  int maxLocalVertex = minVertexIndex + numLocalVertices;
  int vPerProc = numTotalVertices / processors_;
  int locCutsize = 0;
  int totCutsize;
  int numSpanned;
  int arrayLen;
  int totToSend;
  int totToRecv;
  int vertex;
  int endOffset;
  int part;

  int *pVector = &partitionVector[partitionOffsetsVector[pNum]];
  int *array;

  int i;
  int j;
  int ij;

  dynamic_array<int> spannedPart(numParts);
  dynamic_array<int> copyOfReq;

  for (i = 0; i < processors_; ++i)
    send_lens_[i] = 0;

  ds::bit_field sentRequests(numTotalVertices);
  sentRequests.unset();

  // ###
  // compute all the requests for partition
  // vector values of remote vertices
  // ###

  for (i = 0; i < numLocalPins; ++i) {
    ij = localPins[i];

    if (ij < minVertexIndex || ij >= maxLocalVertex) {
      if (!sentRequests(ij)) {
        j = std::min(ij / vPerProc, processors_ - 1);
        data_out_sets_[j]->assign(send_lens_[j]++, ij);
        sentRequests.set(ij);
      }
    }
  }

  // ###
  // compute number of elements to send to other procs
  // ###

  ij = 0;
  for (i = 0; i < processors_; ++i) {
    send_displs_[i] = ij;
    ij += send_lens_[i];
  }

  send_array_.reserve(ij);
  copyOfReq.reserve(ij);

  totToSend = ij;

  ij = 0;
  for (i = 0; i < processors_; ++i) {
    array = data_out_sets_[i]->data();
    arrayLen = send_lens_[i];

    for (j = 0; j < arrayLen; ++j) {
      vertex = array[j];
      send_array_[ij] = vertex;
      copyOfReq[ij++] = vertex;
    }
  }

  MPI_Alltoall(send_lens_.data(), 1, MPI_INT, receive_lens_.data(), 1, MPI_INT,
               comm);

  // ###
  // compute number of elements to receive from other procs
  // ###

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
  // now have received all requests and sent out our requests
  // the reply communication will have the dual dimensions
  // ###

  send_array_.reserve(totToRecv);

  for (i = 0; i < totToRecv; ++i) {
#ifdef DEBUG_HYPERGRAPH
    assert(receive_array_[i] >= minVertexIndex &&
           receive_array_[i] < maxLocalVertex);
#endif

    send_array_[i] = pVector[receive_array_[i] - minVertexIndex];

#ifdef DEBUG_HYPERGRAPH
    assert(send_array_[i] >= 0 && send_array_[i] < numParts);
#endif
  }

  receive_array_.reserve(totToSend);

  MPI_Alltoallv(send_array_.data(), receive_lens_.data(),
                receive_displs_.data(), MPI_INT, receive_array_.data(),
                send_lens_.data(), send_displs_.data(), MPI_INT, comm);

  // ###
  // now the requested vertices are in the copyOfReq data_
  // while their corresponding matchVector values are in
  // the corresponding location in the receive_array_
  // ###

  if (numLocalPins < numTotalVertices / 2) {
    ds::map_from_pos_int<int> storedRequests(numLocalPins);

    for (i = 0; i < totToSend; ++i) {
#ifdef DEBUG_HYPERGRAPH
      assert(receive_array_[i] >= 0 && receive_array_[i] < numParts);
#endif
      storedRequests.insert(copyOfReq[i], receive_array_[i]);
    }

    // ###
    // now procees to compute the contribution to
    // cutsize of the locally held hyperedges
    // ###

    for (i = 0; i < numLocalHedges; ++i) {
      for (j = 0; j < numParts; ++j)
        spannedPart[j] = 0;

      numSpanned = 0;
      endOffset = hEdgeOffsets[i + 1];

      for (j = hEdgeOffsets[i]; j < endOffset; ++j) {
        ij = localPins[j];

        if (ij < minVertexIndex || ij >= maxLocalVertex) {
          part = storedRequests.get(ij);
        } else {
          part = pVector[ij - minVertexIndex];
        }
#ifdef DEBUG_HYPERGRAPH
        assert(part >= 0 && part < numParts);
#endif

        if (spannedPart[part] == 0) {
          spannedPart[part] = 1;
          numSpanned += 1;
        }
#ifdef DEBUG_HYPERGRAPH
        assert(numSpanned >= 1);
#endif
      }

#ifdef DEBUG_HYPERGRAPH
      assert(hEdgeWeights[i] > 0);
      assert(numSpanned > 0);
#endif

      locCutsize += ((numSpanned - 1) * hEdgeWeights[i]);
    }
  } else {
    dynamic_array<int> storedRequests(numTotalVertices);

    for (i = 0; i < totToSend; ++i) {
#ifdef DEBUG_HYPERGRAPH
      assert(receive_array_[i] >= 0 && receive_array_[i] < numParts);
#endif
      storedRequests[copyOfReq[i]] = receive_array_[i];
    }

    // ###
    // now procees to compute the contribution to
    // cutsize of the locally held hyperedges
    // ###

    for (i = 0; i < numLocalHedges; ++i) {
      for (j = 0; j < numParts; ++j)
        spannedPart[j] = 0;

      numSpanned = 0;
      endOffset = hEdgeOffsets[i + 1];

      for (j = hEdgeOffsets[i]; j < endOffset; ++j) {
        ij = localPins[j];

        if (ij < minVertexIndex || ij >= maxLocalVertex) {
          part = storedRequests[ij];
        } else {
          part = pVector[ij - minVertexIndex];
        }
#ifdef DEBUG_HYPERGRAPH
        assert(part >= 0 && part < numParts);
#endif
        if (spannedPart[part] == 0) {
          spannedPart[part] = 1;
          numSpanned += 1;
        }
#ifdef DEBUG_HYPERGRAPH
        assert(numSpanned >= 1);
#endif
      }

#ifdef DEBUG_HYPERGRAPH
      assert(hEdgeWeights[i] > 0);
      assert(numSpanned > 0);
#endif

      locCutsize += ((numSpanned - 1) * hEdgeWeights[i]);
    }
  }

  MPI_Allreduce(&locCutsize, &totCutsize, 1, MPI_INT, MPI_SUM, comm);

  return totCutsize;
}

int ParaHypergraph::checkBalance(int numParts, double balConstraint,
                                 int numPartition, MPI_Comm comm) {
#ifdef DEBUG_HYPERGRAPH
  assert(0 < balConstraint && balConstraint < 0.5);
  assert(numPartition >= 0 && numPartition < numPartitions);
#endif

  int maxPartWt;
  int *pOffset = &partitionVector[partitionOffsetsVector[numPartition]];
  int totWt;

  double avePartWt;

  int i;

  dynamic_array<int> locPWeights(numParts);
  dynamic_array<int> pWeights(numParts);

  for (i = 0; i < numParts; ++i)
    locPWeights[i] = 0;

  for (i = 0; i < numLocalVertices; ++i)
    locPWeights[pOffset[i]] += vWeight[i];

  MPI_Allreduce(locPWeights.data(), pWeights.data(), numParts, MPI_INT,
                MPI_SUM, comm);
  MPI_Allreduce(&localVertexWt, &totWt, 1, MPI_INT, MPI_SUM, comm);

  avePartWt = static_cast<double>(totWt) / numParts;
  maxPartWt = static_cast<int>(floor(avePartWt + avePartWt * balConstraint));

  for (i = 0; i < numParts; ++i)
    if (pWeights[i] > maxPartWt)
      return i;

  return -1;
}

int ParaHypergraph::computeTotalNumPins(MPI_Comm comm) {
  int totalPins;
  MPI_Allreduce(&numLocalPins, &totalPins, 1, MPI_INT, MPI_SUM, comm);
  return totalPins;
}

void ParaHypergraph::checkValidityOfPartitions(int numP) const {
#ifdef DEBUG_HYPERGRAPH
  assert(numPartitions > 0);
#endif

  int i;
  int j = partitionOffsetsVector[numPartitions];

  for (i = 0; i < j; ++i)
    assert(partitionVector[i] >= 0 && partitionVector[i] < numP);
}

void ParaHypergraph::checkPartitions(int numParts, int maxPartWt,
                                     MPI_Comm comm) {
  int i;
  int j;

  int pOffset;

  dynamic_array<int> locPWts(numParts);
  dynamic_array<int> pWts(numParts);

  for (i = 0; i < numPartitions; ++i) {
    for (j = 0; j < numParts; ++j)
      locPWts[j] = 0;

    pOffset = partitionOffsetsVector[i];

    for (j = 0; j < numLocalVertices; ++j) {
      assert(partitionVector[pOffset + j] >= 0 &&
             partitionVector[pOffset + j] < numParts);
      locPWts[partitionVector[pOffset + j]] += vWeight[j];
    }

    MPI_Allreduce(locPWts.data(), pWts.data(), numParts, MPI_INT,
                  MPI_SUM, comm);

    for (j = 0; j < numParts; ++j)
      assert(pWts[j] <= maxPartWt);

    assert(partitionCutsizesVector[i] = calcCutsize(numParts, i, comm));
  }
}

void ParaHypergraph::checkPartitions(int numParts, double constraint,
                                     std::ostream &out, MPI_Comm comm) {
  int i;
  int j;

  int pOffset;
  int cut;
  int totWt;
  int maxWt;
  int maxPartWt;
  double avePartWt;

  char message[512];

  dynamic_array<int> locPWts(numParts);
  dynamic_array<int> pWts(numParts);

  MPI_Allreduce(&localVertexWt, &totWt, 1, MPI_INT, MPI_SUM, comm);

  avePartWt = static_cast<double>(totWt) / numParts;
  maxPartWt = static_cast<int>(floor(avePartWt + avePartWt * constraint));

  for (i = 0; i < numPartitions; ++i) {
    for (j = 0; j < numParts; ++j)
      locPWts[j] = 0;

    pOffset = partitionOffsetsVector[i];

    for (j = 0; j < numLocalVertices; ++j) {
      if (partitionVector[pOffset + j] < 0 ||
          partitionVector[pOffset + j] >= numParts) {
        sprintf(message, "p[%d] - partition vector[%d] = %d\n", rank_,
                minVertexIndex + j, partitionVector[pOffset + j]);
        out << message;
        MPI_Abort(comm, 0);
      }

      locPWts[partitionVector[pOffset + j]] += vWeight[j];
    }

    MPI_Allreduce(locPWts.data(), pWts.data(), numParts, MPI_INT,
                  MPI_SUM, comm);

    if (rank_ == 0) {
      maxWt = 0;
      for (j = 0; j < numParts; ++j)
        if (pWts[j] > maxWt)
          maxWt = pWts[j];

      out << "----- NUM PARTS = " << numParts << std::endl
          << "----- p[" << i << "] largest part weight = " << maxWt << std::endl
          << "----- p[" << i << "] max allowed part weight = " << maxPartWt
          << std::endl;

      if (maxWt <= maxPartWt)
        out << "----- p[" << i << "] satisfies balance constraints" << std::endl;
      else
        out << "----- p[" << i << "] does not satisfy balance constraints"
            << std::endl;
    }

    cut = calcCutsize(numParts, i, comm);

    if (rank_ == 0)
      out << "----- p[" << i << "] k-1 cutsize = " << cut << std::endl;
  }
}

void ParaHypergraph::computeBalanceWarning(int numParts, double constraint,
                                           std::ostream &out, MPI_Comm comm) {
  int i;

  int maxLocVertWt = 0;
  int maxVertWt;
  int maxAllowedVertWt;
  int totWt;

  double avePartWt;

  for (i = 0; i < numLocalVertices; ++i)
    if (vWeight[i] > maxLocVertWt)
      maxLocVertWt = vWeight[i];

  MPI_Reduce(&maxLocVertWt, &maxVertWt, 1, MPI_INT, MPI_MAX, 0, comm);
  MPI_Reduce(&localVertexWt, &totWt, 1, MPI_INT, MPI_SUM, 0, comm);

  if (rank_ == 0) {
    avePartWt = static_cast<double>(totWt) / numParts;
    maxAllowedVertWt = static_cast<int>(floor(avePartWt * constraint));

    if (maxVertWt > maxAllowedVertWt)
      out << "*** Warning! Balance constraint " << constraint
          << " may be too tight ***" << std::endl
          << std::endl;
  }
}

int ParaHypergraph::getNumTotPins(MPI_Comm comm) {
#ifdef DEBUG_HYPERGRAPH
  assert(numLocalPins > 0);
#endif

  int totPins;
  MPI_Allreduce(&numLocalPins, &totPins, 1, MPI_INT, MPI_SUM, comm);
  return totPins;
}

int ParaHypergraph::getNumTotHedges(MPI_Comm comm) {
#ifdef DEBUG_HYPERGRAPH
  assert(numLocalHedges > 0);
#endif

  int totHedges;
  MPI_Allreduce(&numLocalHedges, &totHedges, 1, MPI_INT, MPI_SUM, comm);
  return totHedges;
}

int ParaHypergraph::getExposedHedgeWt(MPI_Comm comm) const {
  int i;
  int ij = 0;
  int totWt;

  for (i = 0; i < numLocalHedges; ++i)
    ij += hEdgeWeights[i];

  MPI_Allreduce(&ij, &totWt, 1, MPI_INT, MPI_SUM, comm);

  return totWt;
}

double ParaHypergraph::getAveVertDeg(MPI_Comm comm) {
  int totPins = getNumTotPins(comm);

  return (static_cast<double>(totPins) / numTotalVertices);
}

double ParaHypergraph::getAveHedgeSize(MPI_Comm comm) {
  int totPins = getNumTotPins(comm);
  int totHedges = getNumTotHedges(comm);

  return (static_cast<double>(totPins) / totHedges);
}

int ParaHypergraph::computeNonConnectedVerts(MPI_Comm comm) {
  ds::bit_field connected(numTotalVertices);
  connected.unset();

  int numTotPins = 0;

  dynamic_array<int> pinsAtProc(processors_);
  dynamic_array<int> receive_displs_(processors_);
  dynamic_array<int> allPins;

  MPI_Gather(&numLocalPins, 1, MPI_INT, pinsAtProc.data(), 1, MPI_INT, 0,
             comm);

  if (rank_ == 0) {
    for (int i = 0; i < processors_; ++i) {
      receive_displs_[i] = numTotPins;
      numTotPins += pinsAtProc[i];
    }

    allPins.reserve(numTotPins);
  }

  assert(localPins.data() != nullptr);
  if (rank_ == 0)
    assert(allPins.data() != nullptr);
  assert(pinsAtProc.data() != nullptr);
  assert(receive_displs_.data() != nullptr);

  MPI_Gatherv(localPins.data(), numLocalPins, MPI_INT, allPins.data(),
              pinsAtProc.data(), receive_displs_.data(), MPI_INT, 0, comm);

  if (rank_ == 0) {
    for (int i = 0; i < numTotPins; ++i)
      connected.set(allPins[i]);

    int numNotConnected = 0;

    for (int i = 0; i < numTotalVertices; ++i)
      if (connected(i) == 0)
        numNotConnected++;

    write_log(rank_, "numConnected = %d", numTotalVertices - numNotConnected);

    return numNotConnected;
  } else
    return 0;
}

#endif
