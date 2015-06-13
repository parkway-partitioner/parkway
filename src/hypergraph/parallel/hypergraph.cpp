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

#include "hypergraph/parallel/hypergraph.hpp"
#include "data_structures/bit_field.hpp"
#include "data_structures/complete_binary_tree.hpp"
#include "data_structures/map_from_pos_int.hpp"
#include "data_structures/new_hyperedge_index_table.hpp"
#include "Log.h"


namespace parkway {
namespace hypergraph {
namespace parallel {

namespace ds = data_structures;

hypergraph::hypergraph(int rank_, int nProcs, int _numLocVerts, int _totVerts,
                       int _minVertIndex, int coarsen, int *wtArray)
    : global_communicator(rank_, nProcs),
      base_hypergraph(_numLocVerts),
      do_not_coarsen(coarsen),
      total_number_of_vertices_(_totVerts),
      minimum_vertex_index_(_minVertIndex),
      vertex_weight_(0) {

  vertex_weights_.set_data(wtArray, number_of_vertices_);
  match_vector_.reserve(number_of_vertices_);
  to_origin_vertex_.reserve(0);

  for (int i = 0; i < number_of_vertices_; ++i) {
    match_vector_[i] = -1;
    vertex_weight_ += vertex_weights_[i];
  }
}

hypergraph::hypergraph(int rank_, int nProcs, int _numLocVerts,
                               int _totVerts, int _minVertIndex, int coarsen,
                               int cut, int *wtArray, int *partArray)
    : global_communicator(rank_, nProcs),
      base_hypergraph(_numLocVerts, 1),
      do_not_coarsen(coarsen),
      total_number_of_vertices_(_totVerts),
      minimum_vertex_index_(_minVertIndex),
      vertex_weight_(0) {

  to_origin_vertex_.reserve(0);
  vertex_weights_.set_data(wtArray, number_of_vertices_);
  partition_vector_.set_data(partArray, number_of_vertices_);
  match_vector_.reserve(number_of_vertices_);
  partition_vector_offsets_.reserve(number_of_partitions_ + 1);
  partition_cuts_.reserve(number_of_partitions_);

  partition_cuts_[0] = cut;
  partition_vector_offsets_[0] = 0;
  partition_vector_offsets_[1] = number_of_vertices_;

  for (int i = 0; i < number_of_vertices_; ++i) {
    match_vector_[i] = -1;
    vertex_weight_ += vertex_weights_[i];
  }
}

hypergraph::hypergraph(int rank_, int nProcs, const char *filename,
                               int dispOption, std::ostream &out, MPI_Comm comm)
    : global_communicator(rank_, nProcs) {
  load_from_file(filename, dispOption, out, comm);
}

hypergraph::hypergraph(int rank_, int nProcs, int numLocVerts,
                               int numLocHedges, int maxHedgeLen,
                               const int *vWeights, const int *hEdgeWts,
                               const int *locPinList, const int *offsets,
                               int dispOption, std::ostream &out, MPI_Comm comm)
    : global_communicator(rank_, nProcs) {

  MPI_Allreduce(&numLocVerts, &total_number_of_vertices_, 1, MPI_INT, MPI_SUM, comm);
  MPI_Scan(&numLocVerts, &minimum_vertex_index_, 1, MPI_INT, MPI_SUM, comm);

  number_of_vertices_ = numLocVerts;
  minimum_vertex_index_ -= number_of_vertices_;
  vertex_weight_ = 0;

  vertex_weights_.reserve(number_of_vertices_);
  to_origin_vertex_.reserve(0);
  match_vector_.reserve(number_of_vertices_);

  for (int i = 0; i < number_of_vertices_; ++i) {
    vertex_weights_[i] = vWeights[i];
    vertex_weight_ += vertex_weights_[i];
    match_vector_[i] = -1;
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
      hyperedge_weights_.assign(hEdgeIndex, hEdgeWts[i]);
      hyperedge_offsets_.assign(hEdgeIndex++, pinCounter);

      for (int j = startOffset; j < endOffset; ++j) {
#ifdef DEBUG_HYPERGRAPH
        assert(locPinList[j] >= 0 && locPinList[j] < numTotalVertices);
#endif
        pin_list_.assign(pinCounter++, locPinList[j]);
      }
    }
  }

  hyperedge_offsets_.assign(hEdgeIndex, pinCounter);

  number_of_pins_ = pinCounter;
  number_of_hyperedges_ = hEdgeIndex;
  do_not_coarsen = 0;
  number_of_partitions_ = 0;

  int i;
  int j;

  if (dispOption > 0) {
    MPI_Reduce(&number_of_pins_, &i, 1, MPI_INT, MPI_SUM, 0, comm);
    MPI_Reduce(&number_of_hyperedges_, &j, 1, MPI_INT, MPI_SUM, 0, comm);

    if (rank_ == 0) {
      out << "|--- Hypergraph (as loaded):" << std::endl;
      out << "| |V| = " << total_number_of_vertices_;
      out << " |E| = " << j;
      out << " |Pins| = " << i << std::endl;
      out << "|" << std::endl;
    }
  }

#ifdef MEM_OPT
  hyperedge_offsets_.reserve(number_of_hyperedges_ + 1);
  hyperedge_weights_.reserve(number_of_hyperedges_);
  pin_list_.reserve(number_of_pins_);
#endif
}

hypergraph::~hypergraph() {}

void hypergraph::load_from_file(const char *filename, int dispOption,
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

  total_number_of_vertices_ = buffer[0];
  number_of_vertices_ = buffer[1];
  hEdgeDataLength = buffer[2];

  minimum_vertex_index_ = (total_number_of_vertices_ / processors_) * rank_;

  to_origin_vertex_.reserve(0);
  vertex_weights_.reserve(number_of_vertices_);
  match_vector_.reserve(number_of_vertices_);
  hEdgeData.reserve(hEdgeDataLength);

  i = sizeof(int) * number_of_vertices_;
  in_stream.read((char *)(vertex_weights_.data()), i);
  if (in_stream.gcount() != i) {
    sprintf(message, "p[%d] could not read in %d vertex elements\n", rank_,
            number_of_vertices_);
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
        << "| |V| = " << total_number_of_vertices_ << " |E| = " << i << std::endl;
  }

  vertex_weight_ = 0;
  for (i = 0; i < number_of_vertices_; ++i) {
    vertex_weight_ += vertex_weights_[i];
    match_vector_[i] = -1;
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
      hyperedge_weights_.assign(hEdgeIndex, hEdgeData[i + 1]);
      hyperedge_offsets_.assign(hEdgeIndex++, pinCounter);

      for (j = 2; j < hEdgeChunk; ++j) {
#ifdef DEBUG_HYPERGRAPH
        assert(hEdgeData[i + j] >= 0 && hEdgeData[i + j] < numTotalVertices);
#endif
        pin_list_.assign(pinCounter++, hEdgeData[i + j]);
      }
    }

    i += hEdgeChunk;
  }

  hyperedge_offsets_.assign(hEdgeIndex, pinCounter);

  number_of_pins_ = pinCounter;
  number_of_hyperedges_ = hEdgeIndex;
  do_not_coarsen = 0;
  number_of_partitions_ = 0;

#ifdef DEBUG_HYPERGRAPH
  for (i = 0; i < numLocalHedges; ++i)
    assert(hEdgeWeights[i] > 0);
#endif

  if (dispOption > 0) {
    MPI_Reduce(&number_of_pins_, &i, 1, MPI_INT, MPI_SUM, 0, comm);
    MPI_Reduce(&number_of_hyperedges_, &j, 1, MPI_INT, MPI_SUM, 0, comm);

    if (rank_ == 0) {
      out << "|--- Hypergraph " << filename << " (as loaded):" << std::endl
          << "| |V| = " << total_number_of_vertices_ << " |E| = " << j
          << " |Pins| = " << i << std::endl
          << "| # Processors = " << processors_ << std::endl
          << "| " << std::endl;
    }
  }

#ifdef MEM_OPT
  hyperedge_offsets_.reserve(number_of_hyperedges_ + 1);
  hyperedge_weights_.reserve(number_of_hyperedges_);
  pin_list_.reserve(number_of_pins_);
#endif
}

void hypergraph::initalize_partition_from_file(const char *filename, int numParts,
                                                        std::ostream &out, MPI_Comm comm) {
  int i;
  int myOffset;
  int numVPerProc = total_number_of_vertices_ / processors_;

  char message[512];

  std::ifstream in_stream;

  number_of_partitions_ = 1;

  partition_vector_.reserve(number_of_vertices_);
  partition_vector_offsets_.reserve(2);
  partition_cuts_.reserve(1);

  partition_vector_offsets_[0] = 0;
  partition_vector_offsets_[1] = number_of_vertices_;

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

  i = number_of_vertices_ * sizeof(int);

  in_stream.read((char *)(partition_vector_.data()), i);
  if (in_stream.gcount() != i) {
    sprintf(message, "p[%d] could not read in %d elements\n", rank_,
            number_of_vertices_);
    out << message;
    MPI_Abort(comm, 0);
  }

  in_stream.close();

  partition_cuts_[0] = calculate_cut_size(numParts, 0, comm);
}

void hypergraph::allocate_hyperedge_memory(int numHedges, int numLocPins) {
#ifdef DEBUG_HYPERGRAPH
  assert(numLocalHedges == numHedges);
  assert(numLocalPins == numLocPins);
#endif

  hyperedge_offsets_.reserve(numHedges + 1);
  hyperedge_weights_.reserve(numHedges);
  pin_list_.reserve(numLocPins);
}

void hypergraph::contract_hyperedges(hypergraph &coarse, MPI_Comm comm) {
  int vFinePerProc;
  int maxLocalVertex = minimum_vertex_index_ + number_of_vertices_;
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

  vFinePerProc = total_number_of_vertices_ / processors_;

  dynamic_array<int> contractedPinList;
  dynamic_array<int> contractedHedgeOffsets;
  dynamic_array<int> contractedHedgeWts;
  dynamic_array<int> origContractedPinList(number_of_pins_);
  dynamic_array<int> copyOfReq;

  for (i = 0; i < processors_; ++i)
    send_lens_[i] = 0;

  ds::bit_field sentRequests(total_number_of_vertices_);
  sentRequests.unset();

  // ###
  // compute all the requests for remote vertex matches
  // ###

  for (i = 0; i < number_of_pins_; ++i) {
    vertex = pin_list_[i];

    if (vertex < minimum_vertex_index_ || vertex >= maxLocalVertex) {
      if (!sentRequests(vertex)) {
        p = std::min(vertex / vFinePerProc, processors_ - 1);
        data_out_sets_[p]->assign(send_lens_[p]++, vertex);
        sentRequests.set(vertex);
      }

      origContractedPinList[i] = -1;
    } else {
      origContractedPinList[i] = match_vector_[vertex - minimum_vertex_index_];
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
    send_array_[i] = match_vector_[receive_array_[i] - minimum_vertex_index_];
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

  if (number_of_pins_ < total_number_of_vertices_ / 2) {
    ds::map_from_pos_int<int> storedRequests(number_of_pins_);

    for (i = 0; i < totalToSend; ++i)
      storedRequests.insert(copyOfReq[i], receive_array_[i]);

    /* contract remaining local pins */

    for (i = 0; i < number_of_pins_; ++i)
      if (origContractedPinList[i] == -1)
        origContractedPinList[i] = storedRequests.get(pin_list_[i]);
  } else {
    dynamic_array<int> nonLocalMatches(total_number_of_vertices_);

    for (i = 0; i < totalToSend; ++i)
      nonLocalMatches[copyOfReq[i]] = receive_array_[i];

    /* contract remaining local pins */

    for (i = 0; i < number_of_pins_; ++i)
      if (origContractedPinList[i] == -1)
        origContractedPinList[i] = nonLocalMatches[pin_list_[i]];
  }

  /* send coarse hyperedges to appropriate processors via hash function */

  for (i = 0; i < number_of_hyperedges_; ++i) {
    j = hyperedge_offsets_[i + 1] - hyperedge_offsets_[i];
    Funct::qsort(0, j - 1, &origContractedPinList[hyperedge_offsets_[i]]);
  }

  contractedPinListLen = 0;
  numContractedHedges = 0;

  for (i = 0; i < number_of_hyperedges_; ++i) {
    contractedHedgeOffsets.assign(numContractedHedges, contractedPinListLen);
    endOffset = hyperedge_offsets_[i + 1];
    startOffset = hyperedge_offsets_[i];

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
      contractedHedgeWts.assign(numContractedHedges++, hyperedge_weights_[i]);

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

  ds::new_hyperedge_index_table table((int)(ceil((double) number_of_hyperedges_ * 1.5)));

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

  coarse.set_number_of_hyperedges(numCoarseHedges);
  coarse.set_number_of_pins(numCoarsePins);
  coarse.allocate_hyperedge_memory(numCoarseHedges, numCoarsePins);

  int *cHedgeOffsets = coarse.hyperedge_offsets();
  int *cPins = coarse.pin_list();
  int *cHedgeWeights = coarse.hyperedge_weights();

  for (i = 0; i < numCoarseHedges; ++i) {
    cHedgeWeights[i] = coarseHedgeWts[i];
    cHedgeOffsets[i] = coarseHedgeOffsets[i];
  }

  cHedgeOffsets[i] = coarseHedgeOffsets[i];

  for (i = 0; i < numCoarsePins; ++i)
    cPins[i] = coarseLocalPins[i];
}

void hypergraph::contractRestrHyperedges(hypergraph &coarse,
                                             MPI_Comm comm) {
  int maxLocalVertex = minimum_vertex_index_ + number_of_vertices_;
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
  dynamic_array<int> origContractedPinList(number_of_pins_);
  dynamic_array<int> copyOfReq;

  MPI_Allgather(&minimum_vertex_index_, 1, MPI_INT, minFineIdxOnProc.data(), 1,
                MPI_INT, comm);

  for (i = 0; i < processors_; ++i)
    procs[i] = i;

  ds::complete_binary_tree <int> vFineToProc(procs.data(),
                                             minFineIdxOnProc.data(), processors_);

  for (i = 0; i < processors_; ++i)
    send_lens_[i] = 0;

  ds::bit_field sentRequests(total_number_of_vertices_);
  sentRequests.unset();

  // ###
  // compute all the requests for remote vertex matches
  // ###

  for (i = 0; i < number_of_pins_; ++i) {
    vertex = pin_list_[i];

    if (vertex < minimum_vertex_index_ || vertex >= maxLocalVertex) {
      if (!sentRequests(vertex)) {
        p = vFineToProc.root_value(vertex);
        data_out_sets_[p]->assign(send_lens_[p]++, vertex);
        sentRequests.set(vertex);
      }

      origContractedPinList[i] = -1;
    } else {
      origContractedPinList[i] = match_vector_[vertex - minimum_vertex_index_];
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
    send_array_[i] = match_vector_[receive_array_[i] - minimum_vertex_index_];
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

  if (number_of_pins_ < total_number_of_vertices_ / 2) {
    ds::map_from_pos_int<int> storedRequests(number_of_pins_);

    for (i = 0; i < totalToSend; ++i)
      storedRequests.insert(copyOfReq[i], receive_array_[i]);

    /* contract remaining local pins */

    for (i = 0; i < number_of_pins_; ++i)
      if (origContractedPinList[i] == -1)
        origContractedPinList[i] = storedRequests.get(pin_list_[i]);
  } else {
    dynamic_array<int> nonLocalMatches(total_number_of_vertices_);

    for (i = 0; i < totalToSend; ++i)
      nonLocalMatches[copyOfReq[i]] = receive_array_[i];

    /* contract remaining local pins */

    for (i = 0; i < number_of_pins_; ++i)
      if (origContractedPinList[i] == -1)
        origContractedPinList[i] = nonLocalMatches[pin_list_[i]];
  }

  /* send coarse hyperedges to appropriate processors via hash function */

  for (i = 0; i < number_of_hyperedges_; ++i) {
    j = hyperedge_offsets_[i + 1] - hyperedge_offsets_[i];
    Funct::qsort(0, j - 1, &origContractedPinList[hyperedge_offsets_[i]]);
  }

  contractedPinListLen = 0;
  numContractedHedges = 0;

  for (i = 0; i < number_of_hyperedges_; ++i) {
    contractedHedgeOffsets.assign(numContractedHedges, contractedPinListLen);
    endOffset = hyperedge_offsets_[i + 1];
    startOffset = hyperedge_offsets_[i];

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
      contractedHedgeWts.assign(numContractedHedges++, hyperedge_weights_[i]);

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

  ds::new_hyperedge_index_table table((int)(ceil((double) number_of_hyperedges_ * 1.1)));

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

  coarse.set_number_of_hyperedges(numCoarseHedges);
  coarse.set_number_of_pins(numCoarsePins);
  coarse.allocate_hyperedge_memory(numCoarseHedges, numCoarsePins);

  int *cHedgeOffsets = coarse.hyperedge_offsets();
  int *cPins = coarse.pin_list();
  int *cHedgeWeights = coarse.hyperedge_weights();

  for (i = 0; i < numCoarseHedges; ++i) {
    cHedgeWeights[i] = coarseHedgeWts[i];
    cHedgeOffsets[i] = coarseHedgeOffsets[i];
  }

  cHedgeOffsets[i] = coarseHedgeOffsets[i];

  for (i = 0; i < numCoarsePins; ++i)
    cPins[i] = coarseLocalPins[i];
}

void hypergraph::project_partitions(hypergraph &coarse, MPI_Comm comm) {
  int totCoarseV = coarse.total_number_of_vertices();
  int numLocCoarseV = coarse.number_of_vertices();
  int minCoarseVindex = coarse.minimum_vertex_index();
  int maxCoarseVindex = minCoarseVindex + numLocCoarseV;
  int numCoarsePartitions = coarse.number_of_partitions();

  int *coarsePartitionVector = coarse.partition_vector();
  int *coarsePartitionOffsets = coarse.partition_offsets();
  int *coarsePartitionCutsizes = coarse.partition_cuts();
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

  number_of_partitions_ = numCoarsePartitions;

  partition_cuts_.reserve(number_of_partitions_);
  partition_vector_offsets_.reserve(number_of_partitions_ + 1);
  partition_vector_.reserve(number_of_partitions_ * number_of_vertices_);

  for (i = 0; i < number_of_partitions_; ++i) {
#ifdef DEBUG_HYPERGRAPH
    assert(coarsePartitionCutsizes[i] > 0);
#endif
    partition_cuts_[i] = coarsePartitionCutsizes[i];
  }

  j = 0;

  for (i = 0; i <= number_of_partitions_; ++i) {
    partition_vector_offsets_[i] = j;
    j += number_of_vertices_;
  }

#ifdef DEBUG_HYPERGRAPH
  for (i = 0; i < partitionOffsetsVector[numPartitions]; ++i)
    partitionVector[i] = -1;
#endif

  for (i = 0; i < number_of_vertices_; ++i) {
#ifdef DEBUG_HYPERGRAPH
    assert(matchVector[i] >= 0 && matchVector[i] < totCoarseV);
#endif

    ij = match_vector_[i];

    if (ij >= minCoarseVindex && ij < maxCoarseVindex) {
      locCoarseVindex = ij - minCoarseVindex;

      for (j = 0; j < number_of_partitions_; ++j) {
        partition_vector_[partition_vector_offsets_[j] + i] =
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

  totToSend = totToRecv * number_of_partitions_;
  send_array_.reserve(totToSend);

  ij = 0;

  for (i = 0; i < totToRecv; ++i) {
#ifdef DEBUG_HYPERGRAPH
    assert(receive_array_[i] >= minCoarseVindex &&
           receive_array_[i] < maxCoarseVindex);
#endif

    vertex = receive_array_[i] - minCoarseVindex;

    for (j = 0; j < number_of_partitions_; ++j) {
      send_array_[ij++] =
          coarsePartitionVector[coarsePartitionOffsets[j] + vertex];
    }
  }

  for (i = 0; i < processors_; ++i) {
    ij = receive_lens_[i];
    receive_lens_[i] = send_lens_[i] * number_of_partitions_;
    send_lens_[i] = ij * number_of_partitions_;
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

    for (j = 0; j < number_of_partitions_; ++j) {
      partition_vector_[partition_vector_offsets_[j] + vertex] = receive_array_[ij++];
    }
  }

#ifdef DEBUG_HYPERGRAPH
  for (i = 0; i < partitionOffsetsVector[numPartitions]; ++i)
    assert(partitionVector[i] != -1);
#endif
}

void hypergraph::reset_vectors() {
  int i;

  for (i = 0; i < number_of_vertices_; ++i)
    match_vector_[i] = -1;

  number_of_partitions_ = 0;

  to_origin_vertex_.reserve(0);
  partition_vector_.reserve(0);
  partition_vector_offsets_.reserve(0);
  partition_cuts_.reserve(0);

  free_memory();
}

void hypergraph::remove_bad_partitions(double cutThreshold) {
#ifdef DEBUG_HYPERGRAPH
  assert(numPartitions > 0);
#endif

  int i;
  int j;

  int bestPartition = 0;
  int bestCut = partition_cuts_[0];
  int acceptedCut;
  int diffInCut;
  int indexIntoOld;
  int indexIntoNew;
  int pSeenBefore;
  int endOffset;
  int numNewPartitions = 0;

  for (i = 1; i < number_of_partitions_; ++i) {
    if (partition_cuts_[i] < bestCut) {
      bestCut = partition_cuts_[i];
      bestPartition = i;
    }
  }

  diffInCut =
      static_cast<int>(floor(static_cast<float>(bestCut) * cutThreshold));
  acceptedCut = bestCut + diffInCut;

  indexIntoNew = 0;
  indexIntoOld = 0;

  for (i = 0; i < number_of_partitions_; ++i) {
    if (partition_cuts_[i] <= acceptedCut) {
      pSeenBefore = 0;

      for (j = 0; j < numNewPartitions; ++j) {
        if (partition_cuts_[j] == partition_cuts_[i])
          pSeenBefore = 1;
      }
      if (pSeenBefore == 0) {
        if (indexIntoOld > indexIntoNew) {
          endOffset = partition_vector_offsets_[i + 1];

          for (; indexIntoOld < endOffset; ++indexIntoOld) {
            partition_vector_[indexIntoNew++] = partition_vector_[indexIntoOld];
          }
        } else {
          indexIntoOld += number_of_vertices_;
          indexIntoNew += number_of_vertices_;
        }

        partition_cuts_[numNewPartitions++] =
            partition_cuts_[i];
      } else {
        indexIntoOld += number_of_vertices_;
      }
    } else {
      indexIntoOld += number_of_vertices_;
    }
  }

  number_of_partitions_ = numNewPartitions;
}

void hypergraph::set_number_of_partitions(int nP) {
  number_of_partitions_ = nP;

  partition_cuts_.reserve(nP);
  partition_vector_offsets_.reserve(nP + 1);

  int j = 0;
  for (int i = 0; i <= number_of_partitions_; ++i) {
    partition_vector_offsets_[i] = j;
    j += number_of_vertices_;
  }

  partition_vector_.reserve(partition_vector_offsets_[number_of_partitions_]);
}

void hypergraph::compute_partition_characteristics(int pNum, int numParts,
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

  int *pVector = &partition_vector_[partition_vector_offsets_[pNum]];

  double avePartWt;

  dynamic_array<int> locPartWeights(numParts);
  dynamic_array<int> partWeights(numParts);

  cut = calculate_cut_size(numParts, pNum, comm);

  for (i = 0; i < numParts; ++i)
    locPartWeights[i] = 0;

  for (i = 0; i < number_of_vertices_; ++i)
    locPartWeights[pVector[i]] += vertex_weights_[i];

  MPI_Reduce(&vertex_weight_, &totHypergraphWeight, 1, MPI_INT, MPI_SUM, 0,
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

void hypergraph::copy_in_partition(const int *partition, int numV,
                                            int nP) {
#ifdef DEBUG_HYPERGRAPH
  assert(numLocalVertices == numV);
  assert(nP >= 0 && nP < numPartitions);
#endif

  int i;

  int endOffset = partition_vector_offsets_[nP + 1];
  int startOffset = partition_vector_offsets_[nP];

  for (i = startOffset; i < endOffset; ++i)
    partition_vector_[i] = partition[i - startOffset];
}

void hypergraph::copy_out_partition(int *partition, int numV,
                                             int nP) const {
#ifdef DEBUG_HYPERGRAPH
  assert(numLocalVertices == numV);
  assert(nP >= 0 && nP < numPartitions);
#endif

  int i;

  int startOffset = partition_vector_offsets_[nP];
  int endOffset = partition_vector_offsets_[nP + 1];

  for (i = startOffset; i < endOffset; ++i)
    partition[i] = partition_vector_[i - startOffset];
}

int hypergraph::keep_best_partition() {
  int i;
  int bestOffset;

  int bestPartition = 0;
  int bestCut = partition_cuts_[0];

  for (i = 1; i < number_of_partitions_; ++i) {
    if (partition_cuts_[i] < bestCut) {
      bestPartition = i;
      bestCut = partition_cuts_[i];
    }
  }

  if (bestPartition != 0) {
    bestOffset = partition_vector_offsets_[bestPartition];

    for (i = 0; i < number_of_vertices_; ++i) {
      partition_vector_[i] = partition_vector_[bestOffset + i];
    }
  }

  partition_cuts_[0] = bestCut;
  number_of_partitions_ = 1;

  return bestCut;
}

void hypergraph::prescribed_vertex_shuffle(int *mapToOrigV, int *prescArray,
                                                    MPI_Comm comm) {
  prescribed_vertex_shuffle(prescArray, number_of_vertices_, comm);
  shift_vertices_to_balance(comm);

  int vPerProc = total_number_of_vertices_ / processors_;
  int maxLocalVertex = minimum_vertex_index_ + number_of_vertices_;
  int totalToRecv;
  int totalToSend;
  int arrayLen;
  int vertex;
  int *array1;
  int *array2;

  int i;
  int j;
  int ij;

  dynamic_array<int> copyOfMapToOrigV(number_of_vertices_);
  dynamic_array<int> askingVertex;

  dynamic_array<dynamic_array<int> *> askingVertices(processors_);

  for (i = 0; i < number_of_vertices_; ++i)
    copyOfMapToOrigV[i] = mapToOrigV[i];

  for (i = 0; i < processors_; ++i) {
    send_lens_[i] = 0;
    askingVertices[i] = new dynamic_array<int>(0);
  }

  /* compute mapToOrigV entries required from un-shuffled hypergraph */

  for (i = 0; i < number_of_vertices_; ++i) {
    vertex = to_origin_vertex_[i];

    if (vertex < minimum_vertex_index_ || vertex >= maxLocalVertex) {
      j = std::min(vertex / vPerProc, processors_ - 1);

      askingVertices[j]->assign(send_lens_[j], i);
      data_out_sets_[j]->assign(send_lens_[j]++, vertex);
    } else {
      mapToOrigV[i] = copyOfMapToOrigV[vertex - minimum_vertex_index_];
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
    send_array_[i] = copyOfMapToOrigV[receive_array_[i] - minimum_vertex_index_];
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

void hypergraph::prescribed_vertex_shuffle(int *prescribedAssignment,
                                                    int nLocVer, MPI_Comm comm) {
#ifdef DEBUG_HYPERGRAPH
  assert(numLocalVertices == nLocVer);
#endif

  int i;

  dynamic_array<int> localVPerProc(processors_);

  for (i = 0; i < processors_; ++i)
    localVPerProc[i] = 0;

  for (i = 0; i < number_of_vertices_; ++i) {
#ifdef DEBUG_HYPERGRAPH
    assert(prescribedAssignment[i] >= 0 && prescribedAssignment[i] < processors_);
#endif
    ++localVPerProc[prescribedAssignment[i]];
  }

  shuffle_vertices(prescribedAssignment, localVPerProc.data(), comm);
}

void hypergraph::shuffle_vertices_by_partition(int nParts, MPI_Comm comm) {
#ifdef DEBUG_HYPERGRAPH
  assert(numPartitions > 0);
  assert(Mod(nParts, processors_) == 0);
#endif

  int i;
  int j;

  int numPartsPerProc = nParts / processors_;

  dynamic_array<int> vToProc(number_of_vertices_);
  dynamic_array<int> localVPerProc(processors_);

  for (i = 0; i < processors_; ++i)
    localVPerProc[i] = 0;

  for (i = 0; i < number_of_vertices_; ++i) {
    j = partition_vector_[i] / numPartsPerProc;
    vToProc[i] = j;
    ++localVPerProc[j];
  }

  shuffle_vertices(vToProc.data(), localVPerProc.data(), comm);
}

void hypergraph::shuffle_vertices_randomly(MPI_Comm comm) {
  int numVerticesEvenlyAllocated;
  int numSpareVertices;
  int totSpareVertices;
  int toProc;

  int i;
  int j;
  int ij;

  dynamic_array<int> vertices(number_of_vertices_);
  dynamic_array<int> vToProc(number_of_vertices_);
  dynamic_array<int> localVPerProc(processors_);
  dynamic_array<int> vSpareToProc;
  dynamic_array<int> indexIntoSpares(processors_);

  for (i = 0; i < processors_; ++i)
    localVPerProc[i] = 0;

  for (i = 0; i < number_of_vertices_; ++i)
    vertices[i] = i;

  Funct::randomPermutation(vertices.data(), number_of_vertices_);

  numSpareVertices = number_of_vertices_ % processors_;
  numVerticesEvenlyAllocated = number_of_vertices_ - numSpareVertices;

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

  for (; i < number_of_vertices_; ++i) {
    j = i - numVerticesEvenlyAllocated;
    toProc = vSpareToProc[indexIntoSpares[rank_] + j];
    vToProc[vertices[i]] = toProc;
    ++localVPerProc[toProc];
  }

  shuffle_vertices(vToProc.data(), localVPerProc.data(), comm);
}

void hypergraph::shuffle_vertices_randomly(int *mapToOrigV, MPI_Comm comm) {
  int numVerticesEvenlyAllocated;
  int numSpareVertices;
  int totSpareVertices;
  int toProc;

  int i;
  int j;
  int ij;

  dynamic_array<int> vertices(number_of_vertices_);
  dynamic_array<int> vToProc(number_of_vertices_);
  dynamic_array<int> localVPerProc(processors_);
  dynamic_array<int> vSpareToProc;
  dynamic_array<int> indexIntoSpares(processors_);

  /* first compute the V->proc map */

  for (i = 0; i < number_of_vertices_; ++i)
    vertices[i] = i;

  Funct::randomPermutation(vertices.data(), number_of_vertices_);

  numSpareVertices = Mod(number_of_vertices_, processors_);
  numVerticesEvenlyAllocated = number_of_vertices_ - numSpareVertices;

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

  j = number_of_vertices_ / processors_;
  for (i = 0; i < processors_; ++i)
    localVPerProc[i] = j;

  for (i = numVerticesEvenlyAllocated; i < number_of_vertices_; ++i) {
    j = i - numVerticesEvenlyAllocated;
    toProc = vSpareToProc[indexIntoSpares[rank_] + j];
    vToProc[vertices[i]] = toProc;
    ++localVPerProc[toProc];
  }

  /* map V->proc is now computed */

  shuffleVerticesAftRandom(vToProc.data(), localVPerProc.data(),
                           mapToOrigV, comm);
}

void hypergraph::shuffle_vertices_randomly(hypergraph &fG, MPI_Comm comm) {
  int numVerticesEvenlyAllocated;
  int numSpareVertices;
  int totSpareVertices;
  int toProc;

  int i;
  int j;
  int ij;

  dynamic_array<int> vertices(number_of_vertices_);
  dynamic_array<int> vToProc(number_of_vertices_);
  dynamic_array<int> localVPerProc(processors_);
  dynamic_array<int> vSpareToProc;
  dynamic_array<int> indexIntoSpares(processors_);

  /* first compute the V->proc map */

  for (i = 0; i < number_of_vertices_; ++i)
    vertices[i] = i;

  Funct::randomPermutation(vertices.data(), number_of_vertices_);

  numSpareVertices = Mod(number_of_vertices_, processors_);
  numVerticesEvenlyAllocated = number_of_vertices_ - numSpareVertices;

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

  j = number_of_vertices_ / processors_;
  for (i = 0; i < processors_; ++i)
    localVPerProc[i] = j;

  for (i = numVerticesEvenlyAllocated; i < number_of_vertices_; ++i) {
    j = i - numVerticesEvenlyAllocated;
    toProc = vSpareToProc[indexIntoSpares[rank_] + j];
    vToProc[vertices[i]] = toProc;
    ++localVPerProc[toProc];
  }

  /* map V->proc is now computed */

  shuffleVerticesAftRandom(vToProc.data(), localVPerProc.data(), fG,
                           comm);
}

void hypergraph::shuffle_vertices(int *vToProc, int *localVPerProc,
                                           MPI_Comm comm) {
#ifdef DEBUG_HYPERGRAPH
  assert(numPartitions == 0 || numPartitions == 1);
#endif

  int i;
  int j;
  int ij;

  int maxLocalVertex = minimum_vertex_index_ + number_of_vertices_;
  int vFinePerProc = total_number_of_vertices_ / processors_;
  int totalToSend;
  int totalToRecv;
  int vertex;
  int startOffset;
  int arrayLen;
  int *array;

  dynamic_array<int> oldIndexToNew(number_of_vertices_);
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

  for (i = 0; i < number_of_vertices_; ++i) {
    oldIndexToNew[i] = newMinVertIndex[vToProc[i]]++;
#ifdef DEBUG_HYPERGRAPH
    assert(oldIndexToNew[i] >= 0 && oldIndexToNew[i] < numTotalVertices);
#endif
  }

  ds::bit_field sentRequests(total_number_of_vertices_);
  sentRequests.unset();

  /* compute vertices required to transform pinlist */

  for (i = 0; i < number_of_pins_; ++i) {
    vertex = pin_list_[i];

    if (vertex < minimum_vertex_index_ || vertex >= maxLocalVertex) {
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
    send_array_[i] = oldIndexToNew[receive_array_[i] - minimum_vertex_index_];
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

  if (number_of_pins_ < total_number_of_vertices_ / 2) {
    ds::map_from_pos_int<int> storedRequests(number_of_pins_);

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

    for (i = 0; i < number_of_pins_; ++i) {
      vertex = pin_list_[i];

      if (vertex >= minimum_vertex_index_ && vertex < maxLocalVertex) {
        pin_list_[i] = oldIndexToNew[vertex - minimum_vertex_index_];
      } else {
        pin_list_[i] = storedRequests.get(vertex);
#ifdef DEBUG_HYPERGRAPH
        assert(localPins[i] >= 0 && localPins[i] < numTotalVertices);
#endif
      }
    }
  } else {
    dynamic_array<int> nonLocalMatches(total_number_of_vertices_);

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

    for (i = 0; i < number_of_pins_; ++i) {
      vertex = pin_list_[i];

      if (vertex >= minimum_vertex_index_ && vertex < maxLocalVertex) {
        pin_list_[i] = oldIndexToNew[vertex - minimum_vertex_index_];
      } else {
        pin_list_[i] = nonLocalMatches[vertex];
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

  if (number_of_partitions_ == 0) {
    totalToSend = Shiftl(number_of_vertices_, 1);

    for (; i < processors_; ++i) {
      send_displs_[i] = j;
      idxIntoSendArray[i] = j;
      send_lens_[i] = Shiftl(localVPerProc[i], 1);
      j += send_lens_[i];
    }
  } else {
    totalToSend = number_of_vertices_ * 3;

    for (; i < processors_; ++i) {
      send_displs_[i] = j;
      idxIntoSendArray[i] = j;
      send_lens_[i] = localVPerProc[i] * 3;
      j += send_lens_[i];
    }
  }

  send_array_.reserve(totalToSend);

  /* compute the send data_ */

  if (number_of_partitions_ == 0) {
    for (i = 0; i < number_of_vertices_; ++i) {
      j = vToProc[i];
      startOffset = idxIntoSendArray[j];
      send_array_[startOffset] = minimum_vertex_index_ + i;
      send_array_[startOffset + 1] = vertex_weights_[i];
      idxIntoSendArray[j] += 2;
    }
  } else {
    for (i = 0; i < number_of_vertices_; ++i) {
      j = vToProc[i];
      startOffset = idxIntoSendArray[j];
      send_array_[startOffset] = minimum_vertex_index_ + i;
      send_array_[startOffset + 1] = vertex_weights_[i];
      send_array_[startOffset + 2] = partition_vector_[i];
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

  number_of_vertices_ = totVperProc[rank_];
  minimum_vertex_index_ = 0;

  for (i = 0; i < rank_; ++i)
    minimum_vertex_index_ += totVperProc[i];

  to_origin_vertex_.reserve(number_of_vertices_);
  vertex_weights_.reserve(number_of_vertices_);
  match_vector_.reserve(number_of_vertices_);

  for (i = 0; i < number_of_vertices_; ++i)
    match_vector_[i] = -1;

  if (number_of_partitions_ > 0) {
    partition_vector_.reserve(number_of_vertices_);
    partition_vector_offsets_[1] = number_of_vertices_;
  }

  i = 0;
  j = 0;
  vertex_weight_ = 0;

  if (number_of_partitions_ == 0) {
    while (i < totalToRecv) {
      to_origin_vertex_[j] = receive_array_[i++];
      vertex_weights_[j++] = receive_array_[i];
      vertex_weight_ += receive_array_[i++];
    }
  } else {
    while (i < totalToRecv) {
      to_origin_vertex_[j] = receive_array_[i++];
      vertex_weights_[j] = receive_array_[i++];
      partition_vector_[j] = receive_array_[i++];
      vertex_weight_ += vertex_weights_[j++];
    }
  }
}

void hypergraph::shuffleVerticesAftRandom(int *vToProc, int *localVPerProc,
                                              int *mapToOrigV, MPI_Comm comm) {
  /* assume that same number of vertices will remain on each processor */

  int i;
  int j;
  int ij;

  int maxLocalVertex = minimum_vertex_index_ + number_of_vertices_;
  int vFinePerProc = total_number_of_vertices_ / processors_;
  int vToOrigVexist = to_origin_vertex_.capacity();
  int totalToSend;
  int totalToRecv;
  int vertex;
  int startOffset;
  int arrayLen;
  int *array;

  dynamic_array<int> oldIndexToNew(number_of_vertices_);
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
  ij = total_number_of_vertices_ / processors_;

  for (i = 0; i < processors_; ++i) {
    minIndexOnProc[i] = j;
    j += ij;
  }

  for (i = 0; i < processors_; ++i)
    minIndexOfMyVerts[i] += minIndexOnProc[i];

  for (i = 0; i < number_of_vertices_; ++i) {
    oldIndexToNew[i] = minIndexOfMyVerts[vToProc[i]]++;
#ifdef DEBUG_HYPERGRAPH
    assert(oldIndexToNew[i] >= 0 && oldIndexToNew[i] < numTotalVertices);
#endif
  }

  ds::bit_field sentRequests(total_number_of_vertices_);
  sentRequests.unset();

  /* compute all the requests for new indices of remote vertices */

  for (i = 0; i < number_of_pins_; ++i) {
    vertex = pin_list_[i];

    if (vertex < minimum_vertex_index_ || vertex >= maxLocalVertex) {
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
    send_array_[i] = oldIndexToNew[receive_array_[i] - minimum_vertex_index_];
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

  if (number_of_pins_ < total_number_of_vertices_ / 2) {
    ds::map_from_pos_int<int> storedRequests(number_of_pins_);

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

    for (i = 0; i < number_of_pins_; ++i) {
      vertex = pin_list_[i];

      if (vertex >= minimum_vertex_index_ && vertex < maxLocalVertex) {
        pin_list_[i] = oldIndexToNew[vertex - minimum_vertex_index_];
      } else {
        pin_list_[i] = storedRequests.get(vertex);
#ifdef DEBUG_HYPERGRAPH
        assert(localPins[i] >= 0 && localPins[i] < numTotalVertices);
#endif
      }
    }
  } else {
    dynamic_array<int> storedRequests(total_number_of_vertices_);

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

    for (i = 0; i < number_of_pins_; ++i) {
      vertex = pin_list_[i];

      if (vertex >= minimum_vertex_index_ && vertex < maxLocalVertex) {
        pin_list_[i] = oldIndexToNew[vertex - minimum_vertex_index_];
      } else {
        pin_list_[i] = storedRequests[vertex];
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
    ij = number_of_partitions_ + 2;
    totalToSend = number_of_vertices_ * ij;
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
    ij = number_of_partitions_ + 3;
    totalToSend = number_of_vertices_ * ij;
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
    for (i = 0; i < number_of_vertices_; ++i) {
      j = vToProc[i];
      startOffset = idxIntoSendArray[j];
      send_array_[startOffset++] = vertex_weights_[i];
      send_array_[startOffset++] = mapToOrigV[i];

      for (ij = 0; ij < number_of_partitions_; ++ij)
        send_array_[startOffset++] =
            partition_vector_[partition_vector_offsets_[ij] + i];

      idxIntoSendArray[j] += (ij + 2);
    }
  } else {
    for (i = 0; i < number_of_vertices_; ++i) {
      j = vToProc[i];
      startOffset = idxIntoSendArray[j];
      send_array_[startOffset++] = to_origin_vertex_[i];
      send_array_[startOffset++] = vertex_weights_[i];
      send_array_[startOffset++] = mapToOrigV[i];

      for (ij = 0; ij < number_of_partitions_; ++ij)
        send_array_[startOffset++] =
            partition_vector_[partition_vector_offsets_[ij] + i];

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
  vertex_weight_ = 0;

  if (!vToOrigVexist) {
    while (i < totalToRecv) {
      vertex_weights_[j] = receive_array_[i];
      vertex_weight_ += receive_array_[i++];
      mapToOrigV[j] = receive_array_[i++];

      for (ij = 0; ij < number_of_partitions_; ++ij)
        partition_vector_[partition_vector_offsets_[ij] + j] = receive_array_[i++];

      ++j;
    }
  } else {
    to_origin_vertex_.reserve(number_of_vertices_);

    while (i < totalToRecv) {
      to_origin_vertex_[j] = receive_array_[i++];
      vertex_weights_[j] = receive_array_[i];
      vertex_weight_ += receive_array_[i++];
      mapToOrigV[j] = receive_array_[i++];

      for (ij = 0; ij < number_of_partitions_; ++ij)
        partition_vector_[partition_vector_offsets_[ij] + j] = receive_array_[i++];

      ++j;
    }
  }
}

void hypergraph::shuffleVerticesAftRandom(int *vToProc, int *localVPerProc,
                                              int *mapToInterV, int *mapToOrigV,
                                              MPI_Comm comm) {
  /* assume that same number of vertices will remain on each processor */

  int i;
  int j;
  int ij;

  int maxLocalVertex = minimum_vertex_index_ + number_of_vertices_;
  int vFinePerProc = total_number_of_vertices_ / processors_;
  int vToOrigVexist = to_origin_vertex_.capacity();
  int totalToSend;
  int totalToRecv;
  int vertex;
  int startOffset;
  int arrayLen;
  int *array;

  dynamic_array<int> oldIndexToNew(number_of_vertices_);
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
  ij = total_number_of_vertices_ / processors_;

  for (i = 0; i < processors_; ++i) {
    minIndexOnProc[i] = j;
    j += ij;
  }

  for (i = 0; i < processors_; ++i)
    minIndexOfMyVerts[i] += minIndexOnProc[i];

  for (i = 0; i < number_of_vertices_; ++i) {
    oldIndexToNew[i] = minIndexOfMyVerts[vToProc[i]]++;
#ifdef DEBUG_HYPERGRAPH
    assert(oldIndexToNew[i] >= 0 && oldIndexToNew[i] < numTotalVertices);
#endif
  }

  /* now need to convert the pinlist and the hyperedge hash keys */

  ds::bit_field sentRequests(total_number_of_vertices_);
  sentRequests.unset();

  /* compute all the requests for new indices of remote vertices  */

  for (i = 0; i < number_of_pins_; ++i) {
    vertex = pin_list_[i];

    if (vertex < minimum_vertex_index_ || vertex >= maxLocalVertex) {
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
    send_array_[i] = oldIndexToNew[receive_array_[i] - minimum_vertex_index_];
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

  if (number_of_pins_ < total_number_of_vertices_ / 2) {
    ds::map_from_pos_int<int> storedRequests(number_of_pins_);

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

    for (i = 0; i < number_of_pins_; ++i) {
      vertex = pin_list_[i];

      if (vertex >= minimum_vertex_index_ && vertex < maxLocalVertex) {
        pin_list_[i] = oldIndexToNew[vertex - minimum_vertex_index_];
      } else {
        pin_list_[i] = storedRequests.get(vertex);
#ifdef DEBUG_HYPERGRAPH
        assert(localPins[i] >= 0 && localPins[i] < numTotalVertices);
#endif
      }
    }
  } else {
    dynamic_array<int> storedRequests(total_number_of_vertices_);

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

    for (i = 0; i < number_of_pins_; ++i) {
      vertex = pin_list_[i];

      if (vertex >= minimum_vertex_index_ && vertex < maxLocalVertex) {
        pin_list_[i] = oldIndexToNew[vertex - minimum_vertex_index_];
      } else {
        pin_list_[i] = storedRequests[vertex];
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
    ij = number_of_partitions_ + 3;
    totalToSend = number_of_vertices_ * ij;

    for (i = 0; i < processors_; ++i) {
      send_displs_[i] = j;
      idxIntoSendArray[i] = j;
      send_lens_[i] = localVPerProc[i] * ij;
      j += send_lens_[i];
    }
  } else {
    ij = number_of_partitions_ + 4;
    totalToSend = number_of_vertices_ * ij;

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
    for (i = 0; i < number_of_vertices_; ++i) {
      j = vToProc[i];
      startOffset = idxIntoSendArray[j];
      send_array_[startOffset++] = vertex_weights_[i];
      send_array_[startOffset++] = mapToInterV[i];
      send_array_[startOffset++] = mapToOrigV[i];

      for (ij = 0; ij < number_of_partitions_; ++ij)
        send_array_[startOffset++] =
            partition_vector_[partition_vector_offsets_[ij] + i];

      idxIntoSendArray[j] += (3 + ij);
    }
  } else {
    for (i = 0; i < number_of_vertices_; ++i) {
      j = vToProc[i];
      startOffset = idxIntoSendArray[j];
      send_array_[startOffset++] = to_origin_vertex_[i];
      send_array_[startOffset++] = vertex_weights_[i];
      send_array_[startOffset++] = mapToInterV[i];
      send_array_[startOffset++] = mapToOrigV[i];

      for (ij = 0; ij < number_of_partitions_; ++ij)
        send_array_[startOffset++] =
            partition_vector_[partition_vector_offsets_[ij] + i];

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
  vertex_weight_ = 0;

  if (!vToOrigVexist) {
    while (i < totalToRecv) {
      vertex_weights_[j] = receive_array_[i];
      vertex_weight_ += receive_array_[i++];
      mapToInterV[j] = receive_array_[i++];
      mapToOrigV[j] = receive_array_[i++];

      for (ij = 0; ij < number_of_partitions_; ++ij)
        partition_vector_[partition_vector_offsets_[ij] + j] = receive_array_[i++];

      ++j;
    }
  } else {
    to_origin_vertex_.reserve(number_of_vertices_);

    while (i < totalToRecv) {
      to_origin_vertex_[j] = receive_array_[i++];
      vertex_weights_[j] = receive_array_[i];
      vertex_weight_ += receive_array_[i++];
      mapToInterV[j] = receive_array_[i++];
      mapToOrigV[j] = receive_array_[i++];

      for (ij = 0; ij < number_of_partitions_; ++ij)
        partition_vector_[partition_vector_offsets_[ij] + j] = receive_array_[i++];

      ++j;
    }
  }
}

void hypergraph::shuffleVerticesAftRandom(int *vToProc, int *localVPerProc,
                                              hypergraph &fineG,
                                              MPI_Comm comm) {
#ifdef DEBUG_HYPERGRAPH
  assert(numPartitions > 0);
#endif

  int i;
  int j;
  int ij;

  int maxLocalVertex = minimum_vertex_index_ + number_of_vertices_;
  int vFinePerProc = total_number_of_vertices_ / processors_;
  int totalToSend;
  int totalToRecv;
  int vertex;
  int startOffset;
  int arrayLen;
  int vToOrigVexist = to_origin_vertex_.capacity();
  int *array;

  dynamic_array<int> oldIndexToNew(number_of_vertices_);
  dynamic_array<int> minIndexOfMyVerts(processors_);
  dynamic_array<int> minIndexOnProc(processors_);
  dynamic_array<int> idxIntoSendArray(processors_);
  dynamic_array<int> copyOfReq;

  /* finer graph structs */

  int numLocalFineVertices = fineG.number_of_vertices();
  int *fineMatchVector = fineG.match_vector();

  for (i = 0; i < processors_; ++i)
    send_lens_[i] = 0;

  MPI_Scan(localVPerProc, minIndexOfMyVerts.data(), processors_, MPI_INT,
           MPI_SUM, comm);

  for (i = 0; i < processors_; ++i)
    minIndexOfMyVerts[i] -= localVPerProc[i];

  j = 0;
  ij = total_number_of_vertices_ / processors_;

  for (i = 0; i < processors_; ++i) {
    minIndexOnProc[i] = j;
    j += ij;
  }

  for (i = 0; i < processors_; ++i)
    minIndexOfMyVerts[i] += minIndexOnProc[i];

  for (i = 0; i < number_of_vertices_; ++i) {
    oldIndexToNew[i] = minIndexOfMyVerts[vToProc[i]]++;
#ifdef DEBUG_HYPERGRAPH
    assert(oldIndexToNew[i] >= 0 && oldIndexToNew[i] < numTotalVertices);
#endif
  }

  /* now need to convert the pinlist and the hyperedge hash keys */

  ds::bit_field sentRequests(total_number_of_vertices_);
  sentRequests.unset();

  /*
     compute all the requests for new indices of remote vertices
     because we need to modify the match vector of the finer graph,
     this needs to include those vertices that are in this match
     vector as well
  */

  for (i = 0; i < number_of_pins_; ++i) {
    vertex = pin_list_[i];

    if (vertex < minimum_vertex_index_ || vertex >= maxLocalVertex) {
      if (!sentRequests(vertex)) {
        j = std::min(vertex / vFinePerProc, processors_ - 1);
        data_out_sets_[j]->assign(send_lens_[j]++, vertex);
        sentRequests.set(vertex);
      }
    }
  }

  for (i = 0; i < numLocalFineVertices; ++i) {
    vertex = fineMatchVector[i];

    if (vertex < minimum_vertex_index_ || vertex >= maxLocalVertex) {
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
    send_array_[i] = oldIndexToNew[receive_array_[i] - minimum_vertex_index_];
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

  if (number_of_pins_ < total_number_of_vertices_ / 2) {
    ds::map_from_pos_int<int> storedRequests(number_of_pins_);

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

    for (i = 0; i < number_of_pins_; ++i) {
      vertex = pin_list_[i];

      if (vertex >= minimum_vertex_index_ && vertex < maxLocalVertex) {
        pin_list_[i] = oldIndexToNew[vertex - minimum_vertex_index_];
      } else {
        pin_list_[i] = storedRequests.get(vertex);
#ifdef DEBUG_HYPERGRAPH
        assert(localPins[i] >= 0 && localPins[i] < numTotalVertices);
#endif
      }
    }

    /* now we will convert the fine graph match vector */

    for (i = 0; i < numLocalFineVertices; ++i) {
      vertex = fineMatchVector[i];

      if (vertex >= minimum_vertex_index_ && vertex < maxLocalVertex) {
        fineMatchVector[i] = oldIndexToNew[vertex - minimum_vertex_index_];
      } else {
        fineMatchVector[i] = storedRequests.get(vertex);
      }
#ifdef DEBUG_HYPERGRAPH
      assert(fineMatchVector[i] >= 0 && fineMatchVector[i] < numTotalVertices);
#endif
    }
  } else {
    dynamic_array<int> storedRequests(total_number_of_vertices_);

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

    for (i = 0; i < number_of_pins_; ++i) {
      vertex = pin_list_[i];

      if (vertex >= minimum_vertex_index_ && vertex < maxLocalVertex) {
        pin_list_[i] = oldIndexToNew[vertex - minimum_vertex_index_];
      } else {
        pin_list_[i] = storedRequests[vertex];
#ifdef DEBUG_HYPERGRAPH
        assert(localPins[i] >= 0 && localPins[i] < numTotalVertices);
#endif
      }
    }

    /* now we will convert the fine graph match vector */

    for (i = 0; i < numLocalFineVertices; ++i) {
      vertex = fineMatchVector[i];

      if (vertex >= minimum_vertex_index_ && vertex < maxLocalVertex) {
        fineMatchVector[i] = oldIndexToNew[vertex - minimum_vertex_index_];
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
    ij = number_of_partitions_ + 2;
    totalToSend = number_of_vertices_ * ij;

    for (i = 0; i < processors_; ++i) {
      send_displs_[i] = j;
      idxIntoSendArray[i] = j;
      send_lens_[i] = localVPerProc[i] * ij;
      j += send_lens_[i];
    }
  } else {
    ij = number_of_partitions_ + 1;
    totalToSend = number_of_vertices_ * ij;

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
    for (i = 0; i < number_of_vertices_; ++i) {
      j = vToProc[i];

      startOffset = idxIntoSendArray[j];
      send_array_[startOffset++] = to_origin_vertex_[i];
      send_array_[startOffset++] = vertex_weights_[i];

      for (ij = 0; ij < number_of_partitions_; ++ij)
        send_array_[startOffset++] =
            partition_vector_[partition_vector_offsets_[ij] + i];

      idxIntoSendArray[j] += (ij + 2);
    }
  } else {
    for (i = 0; i < number_of_vertices_; ++i) {
      j = vToProc[i];

      startOffset = idxIntoSendArray[j];
      send_array_[startOffset++] = vertex_weights_[i];

      for (ij = 0; ij < number_of_partitions_; ++ij)
        send_array_[startOffset++] =
            partition_vector_[partition_vector_offsets_[ij] + i];

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
  vertex_weight_ = 0;

  if (vToOrigVexist > 0) {
    while (i < totalToRecv) {
      to_origin_vertex_[j] = receive_array_[i++];
      vertex_weights_[j] = receive_array_[i];
      vertex_weight_ += receive_array_[i++];

      for (ij = 0; ij < number_of_partitions_; ++ij)
        partition_vector_[partition_vector_offsets_[ij] + j] = receive_array_[i++];

      ++j;
    }
  } else {
    while (i < totalToRecv) {
      vertex_weights_[j] = receive_array_[i];
      vertex_weight_ += receive_array_[i++];

      for (ij = 0; ij < number_of_partitions_; ++ij)
        partition_vector_[partition_vector_offsets_[ij] + j] = receive_array_[i++];

      ++j;
    }
  }
}

void hypergraph::shift_vertices_to_balance(MPI_Comm comm) {
#ifdef DEBUG_HYPERGRAPH
  assert(numPartitions == 0);
#endif

  int i;
  int j;

  int vPerProc = total_number_of_vertices_ / processors_;
  int maxVertexIndex = minimum_vertex_index_ + number_of_vertices_;
  int numMyVertices;

  if (rank_ != processors_ - 1)
    numMyVertices = vPerProc;
  else
    numMyVertices = vPerProc + Mod(total_number_of_vertices_, processors_);

  dynamic_array<int> minNewIndex(processors_);
  dynamic_array<int> maxNewIndex(processors_);

  for (i = 0; i < processors_; ++i) {
    if (i == 0) {
      minNewIndex[i] = 0;
      maxNewIndex[i] = vPerProc;
    } else {
      minNewIndex[i] = maxNewIndex[i - 1];
      if (i == processors_ - 1)
        maxNewIndex[i] = total_number_of_vertices_;
      else
        maxNewIndex[i] = minNewIndex[i] + vPerProc;
    }
  }

  /*
    here distinguish cases where VtoOrigV data_ needs
    to be maintained
  */

  if (to_origin_vertex_.capacity() > 0) {
    j = 0;
    send_array_.reserve(number_of_vertices_ * 3);

    for (i = 0; i < number_of_vertices_; ++i) {
      send_array_[j++] = vertex_weights_[i];
      send_array_[j++] = match_vector_[i];
      send_array_[j++] = to_origin_vertex_[i];
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
          (std::max(number_of_vertices_ - (std::max(maxVertexIndex - maxNewIndex[i], 0) +
                                   std::max(minNewIndex[i] - minimum_vertex_index_, 0)),
               0)) *
          3;
    }
  } else {
    j = 0;
    send_array_.reserve(Shiftl(number_of_vertices_, 1));

    for (i = 0; i < number_of_vertices_; ++i) {
      send_array_[j++] = vertex_weights_[i];
      send_array_[j++] = match_vector_[i];
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
          std::max(number_of_vertices_ - (std::max(maxVertexIndex - maxNewIndex[i], 0) +
                                  std::max(minNewIndex[i] - minimum_vertex_index_, 0)),
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

  number_of_vertices_ = numMyVertices;
  minimum_vertex_index_ = minNewIndex[rank_];

  vertex_weights_.reserve(number_of_vertices_);
  match_vector_.reserve(number_of_vertices_);

  if (to_origin_vertex_.capacity() > 0) {
    to_origin_vertex_.reserve(number_of_vertices_);

    j = 0;
    for (i = 0; i < number_of_vertices_; ++i) {
      vertex_weights_[i] = receive_array_[j++];
      match_vector_[i] = receive_array_[j++];
      to_origin_vertex_[i] = receive_array_[j++];
    }
  } else {
    j = 0;
    for (i = 0; i < number_of_vertices_; ++i) {
      vertex_weights_[i] = receive_array_[j++];
      match_vector_[i] = receive_array_[j++];
    }
  }
}

int hypergraph::calculate_cut_size(int numParts, int pNum, MPI_Comm comm) {
#ifdef DEBUG_HYPERGRAPH
  assert(numPartitions > 0);
  assert(pNum >= 0 && pNum < numPartitions);
#endif

  int maxLocalVertex = minimum_vertex_index_ + number_of_vertices_;
  int vPerProc = total_number_of_vertices_ / processors_;
  int locCutsize = 0;
  int totCutsize;
  int numSpanned;
  int arrayLen;
  int totToSend;
  int totToRecv;
  int vertex;
  int endOffset;
  int part;

  int *pVector = &partition_vector_[partition_vector_offsets_[pNum]];
  int *array;

  int i;
  int j;
  int ij;

  dynamic_array<int> spannedPart(numParts);
  dynamic_array<int> copyOfReq;

  for (i = 0; i < processors_; ++i)
    send_lens_[i] = 0;

  ds::bit_field sentRequests(total_number_of_vertices_);
  sentRequests.unset();

  // ###
  // compute all the requests for partition
  // vector values of remote vertices
  // ###

  for (i = 0; i < number_of_pins_; ++i) {
    ij = pin_list_[i];

    if (ij < minimum_vertex_index_ || ij >= maxLocalVertex) {
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

    send_array_[i] = pVector[receive_array_[i] - minimum_vertex_index_];

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

  if (number_of_pins_ < total_number_of_vertices_ / 2) {
    ds::map_from_pos_int<int> storedRequests(number_of_pins_);

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

    for (i = 0; i < number_of_hyperedges_; ++i) {
      for (j = 0; j < numParts; ++j)
        spannedPart[j] = 0;

      numSpanned = 0;
      endOffset = hyperedge_offsets_[i + 1];

      for (j = hyperedge_offsets_[i]; j < endOffset; ++j) {
        ij = pin_list_[j];

        if (ij < minimum_vertex_index_ || ij >= maxLocalVertex) {
          part = storedRequests.get(ij);
        } else {
          part = pVector[ij - minimum_vertex_index_];
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

      locCutsize += ((numSpanned - 1) * hyperedge_weights_[i]);
    }
  } else {
    dynamic_array<int> storedRequests(total_number_of_vertices_);

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

    for (i = 0; i < number_of_hyperedges_; ++i) {
      for (j = 0; j < numParts; ++j)
        spannedPart[j] = 0;

      numSpanned = 0;
      endOffset = hyperedge_offsets_[i + 1];

      for (j = hyperedge_offsets_[i]; j < endOffset; ++j) {
        ij = pin_list_[j];

        if (ij < minimum_vertex_index_ || ij >= maxLocalVertex) {
          part = storedRequests[ij];
        } else {
          part = pVector[ij - minimum_vertex_index_];
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

      locCutsize += ((numSpanned - 1) * hyperedge_weights_[i]);
    }
  }

  MPI_Allreduce(&locCutsize, &totCutsize, 1, MPI_INT, MPI_SUM, comm);

  return totCutsize;
}

int hypergraph::check_balance(int numParts, double balConstraint,
                                       int numPartition, MPI_Comm comm) {
#ifdef DEBUG_HYPERGRAPH
  assert(0 < balConstraint && balConstraint < 0.5);
  assert(numPartition >= 0 && numPartition < numPartitions);
#endif

  int maxPartWt;
  int *pOffset = &partition_vector_[partition_vector_offsets_[numPartition]];
  int totWt;

  double avePartWt;

  int i;

  dynamic_array<int> locPWeights(numParts);
  dynamic_array<int> pWeights(numParts);

  for (i = 0; i < numParts; ++i)
    locPWeights[i] = 0;

  for (i = 0; i < number_of_vertices_; ++i)
    locPWeights[pOffset[i]] += vertex_weights_[i];

  MPI_Allreduce(locPWeights.data(), pWeights.data(), numParts, MPI_INT,
                MPI_SUM, comm);
  MPI_Allreduce(&vertex_weight_, &totWt, 1, MPI_INT, MPI_SUM, comm);

  avePartWt = static_cast<double>(totWt) / numParts;
  maxPartWt = static_cast<int>(floor(avePartWt + avePartWt * balConstraint));

  for (i = 0; i < numParts; ++i)
    if (pWeights[i] > maxPartWt)
      return i;

  return -1;
}

void hypergraph::check_validity_of_partitions(int numP) const {
#ifdef DEBUG_HYPERGRAPH
  assert(numPartitions > 0);
#endif

  int i;
  int j = partition_vector_offsets_[number_of_partitions_];

  for (i = 0; i < j; ++i)
    assert(partition_vector_[i] >= 0 && partition_vector_[i] < numP);
}

void hypergraph::check_partitions(int numParts, int maxPartWt,
                                           MPI_Comm comm) {
  int i;
  int j;

  int pOffset;

  dynamic_array<int> locPWts(numParts);
  dynamic_array<int> pWts(numParts);

  for (i = 0; i < number_of_partitions_; ++i) {
    for (j = 0; j < numParts; ++j)
      locPWts[j] = 0;

    pOffset = partition_vector_offsets_[i];

    for (j = 0; j < number_of_vertices_; ++j) {
      assert(partition_vector_[pOffset + j] >= 0 &&
             partition_vector_[pOffset + j] < numParts);
      locPWts[partition_vector_[pOffset + j]] += vertex_weights_[j];
    }

    MPI_Allreduce(locPWts.data(), pWts.data(), numParts, MPI_INT,
                  MPI_SUM, comm);

    for (j = 0; j < numParts; ++j)
      assert(pWts[j] <= maxPartWt);

    assert(partition_cuts_[i] = calculate_cut_size(numParts, i, comm));
  }
}

void hypergraph::check_partitions(int numParts, double constraint,
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

  MPI_Allreduce(&vertex_weight_, &totWt, 1, MPI_INT, MPI_SUM, comm);

  avePartWt = static_cast<double>(totWt) / numParts;
  maxPartWt = static_cast<int>(floor(avePartWt + avePartWt * constraint));

  for (i = 0; i < number_of_partitions_; ++i) {
    for (j = 0; j < numParts; ++j)
      locPWts[j] = 0;

    pOffset = partition_vector_offsets_[i];

    for (j = 0; j < number_of_vertices_; ++j) {
      if (partition_vector_[pOffset + j] < 0 ||
          partition_vector_[pOffset + j] >= numParts) {
        sprintf(message, "p[%d] - partition vector[%d] = %d\n", rank_,
                minimum_vertex_index_ + j, partition_vector_[pOffset + j]);
        out << message;
        MPI_Abort(comm, 0);
      }

      locPWts[partition_vector_[pOffset + j]] += vertex_weights_[j];
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

    cut = calculate_cut_size(numParts, i, comm);

    if (rank_ == 0)
      out << "----- p[" << i << "] k-1 cutsize = " << cut << std::endl;
  }
}

void hypergraph::compute_balance_warnings(int numParts, double constraint,
                                                   std::ostream &out, MPI_Comm comm) {
  int i;

  int maxLocVertWt = 0;
  int maxVertWt;
  int maxAllowedVertWt;
  int totWt;

  double avePartWt;

  for (i = 0; i < number_of_vertices_; ++i)
    if (vertex_weights_[i] > maxLocVertWt)
      maxLocVertWt = vertex_weights_[i];

  MPI_Reduce(&maxLocVertWt, &maxVertWt, 1, MPI_INT, MPI_MAX, 0, comm);
  MPI_Reduce(&vertex_weight_, &totWt, 1, MPI_INT, MPI_SUM, 0, comm);

  if (rank_ == 0) {
    avePartWt = static_cast<double>(totWt) / numParts;
    maxAllowedVertWt = static_cast<int>(floor(avePartWt * constraint));

    if (maxVertWt > maxAllowedVertWt)
      out << "*** Warning! Balance constraint " << constraint
          << " may be too tight ***" << std::endl
          << std::endl;
  }
}

int hypergraph::total_number_of_pins(MPI_Comm comm) {
#ifdef DEBUG_HYPERGRAPH
  assert(numLocalPins > 0);
#endif

  int totPins;
  MPI_Allreduce(&number_of_pins_, &totPins, 1, MPI_INT, MPI_SUM, comm);
  return totPins;
}

int hypergraph::total_number_of_hyperedges(MPI_Comm comm) {
#ifdef DEBUG_HYPERGRAPH
  assert(numLocalHedges > 0);
#endif

  int totHedges;
  MPI_Allreduce(&number_of_hyperedges_, &totHedges, 1, MPI_INT, MPI_SUM, comm);
  return totHedges;
}

int hypergraph::exposed_hyperedge_weight(MPI_Comm comm) const {
  int i;
  int ij = 0;
  int totWt;

  for (i = 0; i < number_of_hyperedges_; ++i)
    ij += hyperedge_weights_[i];

  MPI_Allreduce(&ij, &totWt, 1, MPI_INT, MPI_SUM, comm);

  return totWt;
}

double hypergraph::average_vertex_degree(MPI_Comm comm) {
  int totPins = total_number_of_pins(comm);

  return (static_cast<double>(totPins) / total_number_of_vertices_);
}

double hypergraph::average_hyperedge_size(MPI_Comm comm) {
  int totPins = total_number_of_pins(comm);
  int totHedges = total_number_of_hyperedges(comm);

  return (static_cast<double>(totPins) / totHedges);
}

}  // namespace parallel
}  // namespace hypergraph
}  // namespace parkway
