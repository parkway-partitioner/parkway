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
#include "utility/sorting.hpp"
#include "Log.h"

#include <bitset>

namespace parkway {
namespace parallel {
namespace ds = parkway::data_structures;
hypergraph::hypergraph(int rank, int number_of_processors,
                       int number_of_local_vertices, int total_vertices,
                       int minimum_vertex_index, int coarsen,
                       ds::dynamic_array<int> weights,
                       int display_option)
    : global_communicator(rank, number_of_processors, display_option),
      base_hypergraph(number_of_local_vertices),
      do_not_coarsen(coarsen),
      total_number_of_vertices_(total_vertices),
      minimum_vertex_index_(minimum_vertex_index),
      vertex_weight_(0) {
  vertex_weights_ = weights;
  match_vector_.assign(number_of_vertices_, -1);
  for (int i = 0; i < number_of_vertices_; ++i) {
    vertex_weight_ += vertex_weights_[i];
  }
}


hypergraph::hypergraph(int rank, int number_of_processors,
                       int number_of_local_vertices, int total_vertices,
                       int minimum_vertex_index, int coarsen, int cut,
                       ds::dynamic_array<int> weight,
                       ds::dynamic_array<int> part_array,
                       int display_option)
    : global_communicator(rank, number_of_processors, display_option),
      base_hypergraph(number_of_local_vertices, 1),
      do_not_coarsen(coarsen),
      total_number_of_vertices_(total_vertices),
      minimum_vertex_index_(minimum_vertex_index),
      vertex_weight_(0) {
  vertex_weights_ = weight;
  partition_vector_ = part_array;
  match_vector_.assign(number_of_vertices_, -1);
  partition_vector_offsets_.reserve(number_of_partitions_ + 1);
  partition_cuts_.reserve(number_of_partitions_);

  partition_cuts_[0] = cut;
  partition_vector_offsets_[0] = 0;
  partition_vector_offsets_[1] = number_of_vertices_;

  for (const int &weight : vertex_weights_) {
    vertex_weight_ += weight;
  }
}

hypergraph::hypergraph(int rank, int number_of_processors, const char *filename,
                       int dispOption, std::ostream &out, MPI_Comm comm)
    : global_communicator(rank, number_of_processors, dispOption) {
  load_from_file(filename, out, comm);
}

hypergraph::hypergraph(int rank, int number_of_processors,
                       int number_of_local_vertices, int number_of_local_hedges,
                       int max_hyperedge_length,
                       ds::dynamic_array<int> vertex_weights,
                       ds::dynamic_array<int> hyperedge_weights,
                       ds::dynamic_array<int> loc_pin_list,
                       ds::dynamic_array<int> hyperedge_offsets,
                       int display_option, std::ostream &out, MPI_Comm comm)
    : global_communicator(rank_, number_of_processors, display_option) {
  MPI_Allreduce(&number_of_local_vertices, &total_number_of_vertices_, 1,
                MPI_INT, MPI_SUM, comm);
  MPI_Scan(&number_of_local_vertices, &minimum_vertex_index_, 1, MPI_INT,
           MPI_SUM, comm);

  number_of_vertices_ = number_of_local_vertices;
  minimum_vertex_index_ -= number_of_vertices_;
  vertex_weight_ = 0;

  vertex_weights_.reserve(number_of_vertices_);
  match_vector_.reserve(number_of_vertices_);

  for (int i = 0; i < number_of_vertices_; ++i) {
    vertex_weights_[i] = vertex_weights[i];
    vertex_weight_ += vertex_weights_[i];
    match_vector_[i] = -1;
  }

  /* hyperedges stored in pinlist */
  int hyperedge_index = 0;
  int pin_counter = 0;
  for (int i = 0; i < number_of_local_hedges; ++i) {
    int start_offset = hyperedge_offsets[i];
    int end_offset = hyperedge_offsets[i + 1];
    int length = end_offset - start_offset;

    if (length > 1) {
      hyperedge_weights_[hyperedge_index] = hyperedge_weights[i];
      hyperedge_offsets_[hyperedge_index++] = pin_counter;
      for (int j = start_offset; j < end_offset; ++j) {
        pin_list_[pin_counter++] = loc_pin_list[j];
      }
    }
  }

  hyperedge_offsets_[hyperedge_index] = pin_counter;

  number_of_pins_ = pin_counter;
  number_of_hyperedges_ = hyperedge_index;
  do_not_coarsen = 0;
  number_of_partitions_ = 0;

  if (display_option > 0) {
    int i;
    int j;
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
}

hypergraph::~hypergraph() {
}

void hypergraph::load_from_file(const char *filename, std::ostream &out,
                                MPI_Comm comm) {
  // Format filename so that each processor loads the correct part.
  char my_file[512];
  sprintf(my_file, "%s-%d", filename, rank_);
  std::ifstream in_stream(my_file, std::ifstream::in | std::ifstream::binary);

  // Print message if the file could not be opened.
  char message[512];
  if (!in_stream.is_open()) {
    out << "p[" << rank_ << "] could not open " << my_file << std::endl;
    MPI_Abort(comm, 0);
  }

  // Read the hypergraph metadata -- total number of vertices, vertices in this
  // part and number of hyperedges.
  int buffer[3];
  int buffer_size = sizeof(int) * 3;
  in_stream.read((char *)(&buffer[0]), buffer_size);
  if (in_stream.gcount() != buffer_size) {
    out << "p[" << rank_ << "] could not read in metadata (three integers)"
        << std::endl;
    in_stream.close();
    MPI_Abort(comm, 0);
  }

  total_number_of_vertices_ = buffer[0];
  number_of_vertices_ = buffer[1];
  int hyperedge_data_length = buffer[2];

  minimum_vertex_index_ = (total_number_of_vertices_ / processors_) * rank_;

  // Read in vertex weights and hyperedge data
  if (!vertex_weights_.read_from(in_stream, number_of_vertices_)) {
    out << "p[" << rank_ << "] could not read in " << number_of_vertices_
        << " vertex elements" << std::endl;
    in_stream.close();
    MPI_Abort(comm, 0);
  }

  dynamic_array<int> hyperedge_data(hyperedge_data_length);
  if (!hyperedge_data.read_from(in_stream, hyperedge_data_length)) {
    out << "p[" << rank_ << "] could not read in " << hyperedge_data_length
        << " hyperedge elements" << std::endl;
    in_stream.close();
    MPI_Abort(comm, 0);
  }
  in_stream.close();


  check_vertex_and_hyperedge_lengths(hyperedge_data_length, hyperedge_data,
                                     comm);

  match_vector_.assign(number_of_vertices_, -1);
  vertex_weight_ = 0;
  for (const auto &weight : vertex_weights_) {
    vertex_weight_ += weight;
  }

  load_data_from_blocks(hyperedge_data_length, hyperedge_data);
  check_loaded_vertex_and_hyperedge_lengths();
}

void hypergraph::initalize_partition_from_file(const char *filename,
                                               int numParts, std::ostream &out,
                                               MPI_Comm comm) {
  int i;
  int vertices_per_processor = total_number_of_vertices_ / processors_;

  char message[512];

  number_of_partitions_ = 1;

  partition_vector_.resize(number_of_vertices_);
  partition_cuts_.resize(1);

  partition_vector_offsets_.resize(2);
  partition_vector_offsets_[0] = 0;
  partition_vector_offsets_[1] = number_of_vertices_;


  std::ifstream in_stream(filename, std::ios::in | std::ios::binary);
  if (!in_stream.is_open()) {
    out << "p[" << rank_ << "] could not open partition file '" << filename
        << "'" << std::endl;
    MPI_Abort(comm, 0);
  }

  int offset_on_this_rank = rank_ * vertices_per_processor;
  in_stream.seekg(offset_on_this_rank * sizeof(int), std::ifstream::beg);

  if (!partition_vector_.read_from(in_stream, number_of_vertices_)) {
    out << "p[" << rank_ << "] could not read in " << number_of_vertices_
        << " elements" << std::endl;
    MPI_Abort(comm, 0);
  }

  in_stream.close();
  partition_cuts_[0] = calculate_cut_size(numParts, 0, comm);
}

void hypergraph::allocate_hyperedge_memory(int numHedges, int numLocPins) {
  hyperedge_offsets_.reserve(numHedges + 1);
  hyperedge_weights_.reserve(numHedges);
  pin_list_.reserve(numLocPins);
}

void hypergraph::contract_hyperedges(hypergraph &coarse, MPI_Comm comm) {
  int max_local_vertex = minimum_vertex_index_ + number_of_vertices_;
  int vertices_per_proc_fine = total_number_of_vertices_ / processors_;

  dynamic_array<int> original_contracted_pin_list(number_of_pins_, -1);
  compute_requests_for_remote_vertex_matches(vertices_per_proc_fine);

  // compute number of elements to send to other processors
  ds::dynamic_array<int> copy_of_requests;
  int total_to_send = compute_number_of_elements_to_send(copy_of_requests);
  MPI_Alltoall(send_lens_.data(), 1, MPI_INT, receive_lens_.data(), 1, MPI_INT,
               comm);

  // compute number of elements to receive from other processors
  int total_to_receive = compute_number_of_elements_to_receive();
  MPI_Alltoallv(send_array_.data(), send_lens_.data(),
                send_displs_.data(), MPI_INT, receive_array_.data(),
                receive_lens_.data(), receive_displs_.data(), MPI_INT, comm);

  // now have received all requests and sent out our requests the reply
  // communication will have the dual dimensions
  send_array_.resize(total_to_receive);
  for (int i = 0; i < total_to_receive; ++i) {
    send_array_[i] = match_vector_[receive_array_[i] - minimum_vertex_index_];
  }

  receive_array_.resize(total_to_send);
  MPI_Alltoallv(send_array_.data(), receive_lens_.data(),
                receive_displs_.data(), MPI_INT, receive_array_.data(),
                send_lens_.data(), send_displs_.data(), MPI_INT, comm);

  // Requested vertices are in the copy_of_requests while their corresponding
  // match_vector values are in the corresponding location in the receive_array_
  //
  // Depending on number of local pins, choose format for storing non-local
  // vertices
  choose_non_local_vertices_format(original_contracted_pin_list);

  // send coarse hyperedges to appropriate processors via hash function
  send_coarse_hyperedges(original_contracted_pin_list, total_to_send,
                         total_to_receive);

  // Should have received all hyperedges destined for processor now build the
  // coarse hypergraph pin-list
  int table_size = static_cast<int>(
      ceil(static_cast<double>(number_of_hyperedges_) * 1.5));
  ds::new_hyperedge_index_table table(table_size);

  // run through the recv data_
  // for each encountered hyperedge:
  // - compute hash-key
  // - check for duplicates in the hash table
  // - if duplicate not found, add to list
  // - if duplicate found, increment its weight
  process_new_hyperedges(coarse, table);
}


void hypergraph::contractRestrHyperedges(hypergraph &coarse, MPI_Comm comm) {
  int max_local_vertex = minimum_vertex_index_ + number_of_vertices_;

  dynamic_array<int> min_fine_index_on_processor(processors_);
  dynamic_array<int> original_contracted_pin_list(number_of_pins_);

  MPI_Allgather(&minimum_vertex_index_, 1, MPI_INT,
                min_fine_index_on_processor.data(), 1, MPI_INT, comm);

  dynamic_array<int> processors(processors_);
  for (int i = 0; i < processors_; ++i) {
    processors[i] = i;
  }

  ds::complete_binary_tree <int> fine_vertex_to_processor(
      processors.data(), min_fine_index_on_processor.data(), processors_);

  compute_requests_for_remote_vertex_matches(0, &fine_vertex_to_processor);

  // Compute number of elements to send to other processors
  ds::dynamic_array<int> copy_of_requests;
  int total_to_send = compute_number_of_elements_to_send(copy_of_requests);
  MPI_Alltoall(send_lens_.data(), 1, MPI_INT, receive_lens_.data(), 1, MPI_INT,
               comm);

  // Compute number of elements to receive from other processors
  int total_to_receive = compute_number_of_elements_to_receive();
  MPI_Alltoallv(send_array_.data(), send_lens_.data(),
                send_displs_.data(), MPI_INT, receive_array_.data(),
                receive_lens_.data(), receive_displs_.data(), MPI_INT, comm);

  // now have received all requests and sent out our requests the reply
  // communication will have the dual dimensions
  send_array_.resize(total_to_receive);
  for (int i = 0; i < total_to_receive; ++i) {
    send_array_[i] = match_vector_[receive_array_[i] - minimum_vertex_index_];
  }

  receive_array_.resize(total_to_send);
  MPI_Alltoallv(send_array_.data(), receive_lens_.data(),
                receive_displs_.data(), MPI_INT, receive_array_.data(),
                send_lens_.data(), send_displs_.data(), MPI_INT, comm);

  // Requested vertices are in the copy_of_requests data while their
  // corresponding matchVector values are in the corresponding location in the
  // receieve_array_
  choose_non_local_vertices_format(original_contracted_pin_list);

  send_coarse_hyperedges(original_contracted_pin_list, total_to_send,
                         total_to_receive);

  // Should have received all hyperedges destined for processor now build the
  // coarse hypergraph pin-list
  int table_size = static_cast<int>(
      ceil(static_cast<double>(number_of_hyperedges_) * 1.1));
  ds::new_hyperedge_index_table table(table_size);

  // run through the recv data_
  // for each encountered hyperedge:
  // - compute hash-key
  // - check for duplicates in the hash table
  // - if duplicate not found, add to list
  // - if duplicate found, increment its weight

  // NEED TO MODIFY THE RESTR HEDGE CONTRACTION TO CORRESPOND WITH NORMAL.
  process_new_hyperedges(coarse, table);
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
  int end_offset;
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
          end_offset = partition_vector_offsets_[i + 1];

          for (; indexIntoOld < end_offset; ++indexIntoOld) {
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

  int end_offset = partition_vector_offsets_[nP + 1];
  int start_offset = partition_vector_offsets_[nP];

  for (i = start_offset; i < end_offset; ++i)
    partition_vector_[i] = partition[i - start_offset];
}

void hypergraph::copy_out_partition(int *partition, int numV,
                                             int nP) const {
#ifdef DEBUG_HYPERGRAPH
  assert(numLocalVertices == numV);
  assert(nP >= 0 && nP < numPartitions);
#endif

  int i;

  int start_offset = partition_vector_offsets_[nP];
  int end_offset = partition_vector_offsets_[nP + 1];

  for (i = start_offset; i < end_offset; ++i)
    partition[i] = partition_vector_[i - start_offset];
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
  int max_local_vertex = minimum_vertex_index_ + number_of_vertices_;
  int total_to_receive;
  int total_to_send;
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

    if (vertex < minimum_vertex_index_ || vertex >= max_local_vertex) {
      j = std::min(vertex / vPerProc, processors_ - 1);

      askingVertices[j]->assign(send_lens_[j], i);
      data_out_sets_[j]->assign(send_lens_[j]++, vertex);
    } else {
      mapToOrigV[i] = copyOfMapToOrigV[vertex - minimum_vertex_index_];
    }
  }

  /* compute number of elements to send to other processors */

  j = 0;
  for (i = 0; i < processors_; ++i) {
    send_displs_[i] = j;
    j += send_lens_[i];
  }

  send_array_.reserve(j);
  askingVertex.reserve(j);
  total_to_send = j;

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
  assert(j == total_to_send);
#endif

  MPI_Alltoall(send_lens_.data(), 1, MPI_INT, receive_lens_.data(), 1, MPI_INT,
               comm);

  /* compute number of elements to receive from other processors */

  j = 0;
  for (i = 0; i < processors_; ++i) {
    receive_displs_[i] = j;
    j += receive_lens_[i];
  }

  receive_array_.reserve(j);
  total_to_receive = j;

  MPI_Alltoallv(send_array_.data(), send_lens_.data(),
                send_displs_.data(), MPI_INT, receive_array_.data(),
                receive_lens_.data(), receive_displs_.data(), MPI_INT, comm);

  /*
    now have received all requests and sent out our requests
    the reply communication will have the dual dimensions
  */

  send_array_.reserve(total_to_receive);

  for (i = 0; i < total_to_receive; ++i) {
#ifdef DEBUG_HYPERGRAPH
    assert(receive_array_[i] >= minVertexIndex &&
           receive_array_[i] < max_local_vertex);
#endif
    send_array_[i] = copyOfMapToOrigV[receive_array_[i] - minimum_vertex_index_];
#ifdef DEBUG_HYPERGRAPH
    assert(send_array_[i] >= 0 && send_array_[i] < numTotalVertices);
#endif
  }

  receive_array_.reserve(total_to_send);

  MPI_Alltoallv(send_array_.data(), receive_lens_.data(),
                receive_displs_.data(), MPI_INT, receive_array_.data(),
                send_lens_.data(), send_displs_.data(), MPI_INT, comm);

  /*
     now the requested vertices are in the copy_of_requests data_
     while their corresponding matchVector values are in
     the corresponding location in the receive_array_
  */

  for (i = 0; i < total_to_send; ++i) {
#ifdef DEBUG_HYPERGRAPH
    assert(receive_array_[i] >= 0 && receive_array_[i] < numTotalVertices);
#endif
    mapToOrigV[askingVertex[i]] = receive_array_[i];
  }

  for (i = 0; i < processors_; ++i)
    dynamic_memory::delete_pointer<dynamic_array<int> >(askingVertices[i]);
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

  int max_local_vertex = minimum_vertex_index_ + number_of_vertices_;
  int vertices_per_proc_fine = total_number_of_vertices_ / processors_;
  int total_to_send;
  int total_to_receive;
  int vertex;
  int start_offset;
  int arrayLen;
  int *array;

  dynamic_array<int> oldIndexToNew(number_of_vertices_);
  dynamic_array<int> newMinVertIndex(processors_);
  dynamic_array<int> totVperProc(processors_);
  dynamic_array<int> minNewIndexOnProc(processors_);

  dynamic_array<int> idxIntoSendArray(processors_);
  dynamic_array<int> copy_of_requests;

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

  ds::bit_field sent_requests(total_number_of_vertices_);
  sent_requests.unset();

  /* compute vertices required to transform pinlist */

  for (i = 0; i < number_of_pins_; ++i) {
    vertex = pin_list_[i];

    if (vertex < minimum_vertex_index_ || vertex >= max_local_vertex) {
      if (!sent_requests(vertex)) {
        j = std::min(vertex / vertices_per_proc_fine, processors_ - 1);
        data_out_sets_[j]->assign(send_lens_[j]++, vertex);
        sent_requests.set(vertex);
      }
    }
  }

  /* compute number of elements to send to other processors */

  j = 0;
  for (i = 0; i < processors_; ++i) {
    send_displs_[i] = j;
    j += send_lens_[i];
  }

  send_array_.reserve(j);
  copy_of_requests.reserve(j);
  total_to_send = j;

  j = 0;
  for (ij = 0; ij < processors_; ++ij) {
    array = data_out_sets_[ij]->data();
    arrayLen = send_lens_[ij];

    for (i = 0; i < arrayLen; ++i) {
      vertex = array[i];
      send_array_[j] = vertex;
      copy_of_requests[j++] = vertex;
    }
  }

#ifdef DEBUG_HYPERGRAPH
  assert(j == total_to_send);
#endif

  MPI_Alltoall(send_lens_.data(), 1, MPI_INT, receive_lens_.data(), 1, MPI_INT,
               comm);

  /* compute number of elements to receive from other processors */

  j = 0;
  for (i = 0; i < processors_; ++i) {
    receive_displs_[i] = j;
    j += receive_lens_[i];
  }

  receive_array_.reserve(j);
  total_to_receive = j;

  MPI_Alltoallv(send_array_.data(), send_lens_.data(),
                send_displs_.data(), MPI_INT, receive_array_.data(),
                receive_lens_.data(), receive_displs_.data(), MPI_INT, comm);

  /*
    now have received all requests and sent out our requests
    the reply communication will have the dual dimensions
  */

  send_array_.reserve(total_to_receive);

  for (i = 0; i < total_to_receive; ++i) {
#ifdef DEBUG_HYPERGRAPH
    assert(receive_array_[i] >= minVertexIndex &&
           receive_array_[i] < max_local_vertex);
#endif
    send_array_[i] = oldIndexToNew[receive_array_[i] - minimum_vertex_index_];
#ifdef DEBUG_HYPERGRAPH
    assert(send_array_[i] >= 0 && send_array_[i] < numTotalVertices);
#endif
  }

  receive_array_.reserve(total_to_send);

  MPI_Alltoallv(send_array_.data(), receive_lens_.data(),
                receive_displs_.data(), MPI_INT, receive_array_.data(),
                send_lens_.data(), send_displs_.data(), MPI_INT, comm);

  /*
   now the requested vertices are in the copy_of_requests data_
   while their corresponding matchVector values are in
   the corresponding location in the receive_array_
  */

  if (number_of_pins_ < total_number_of_vertices_ / 2) {
    ds::map_from_pos_int<int> stored_requests(number_of_pins_);

    for (i = 0; i < total_to_send; ++i) {
#ifdef DEBUG_HYPERGRAPH
      assert(receive_array_[i] >= 0 && receive_array_[i] < numTotalVertices);
#endif
      stored_requests.insert(copy_of_requests[i], receive_array_[i]);
    }

#ifdef DEBUG_HYPERGRAPH
    for (i = 0; i < numLocalVertices; ++i)
      assert(oldIndexToNew[i] >= 0 && oldIndexToNew[i] < numTotalVertices);
#endif

    for (i = 0; i < number_of_pins_; ++i) {
      vertex = pin_list_[i];

      if (vertex >= minimum_vertex_index_ && vertex < max_local_vertex) {
        pin_list_[i] = oldIndexToNew[vertex - minimum_vertex_index_];
      } else {
        pin_list_[i] = stored_requests.get(vertex);
#ifdef DEBUG_HYPERGRAPH
        assert(localPins[i] >= 0 && localPins[i] < numTotalVertices);
#endif
      }
    }
  } else {
    dynamic_array<int> non_local_matches(total_number_of_vertices_);

    for (i = 0; i < total_to_send; ++i) {
#ifdef DEBUG_HYPERGRAPH
      assert(receive_array_[i] >= 0 && receive_array_[i] < numTotalVertices);
#endif
      non_local_matches[copy_of_requests[i]] = receive_array_[i];
    }

#ifdef DEBUG_HYPERGRAPH
    for (i = 0; i < numLocalVertices; ++i)
      assert(oldIndexToNew[i] >= 0 && oldIndexToNew[i] < numTotalVertices);
#endif

    for (i = 0; i < number_of_pins_; ++i) {
      vertex = pin_list_[i];

      if (vertex >= minimum_vertex_index_ && vertex < max_local_vertex) {
        pin_list_[i] = oldIndexToNew[vertex - minimum_vertex_index_];
      } else {
        pin_list_[i] = non_local_matches[vertex];
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
    total_to_send = Shiftl(number_of_vertices_, 1);

    for (; i < processors_; ++i) {
      send_displs_[i] = j;
      idxIntoSendArray[i] = j;
      send_lens_[i] = Shiftl(localVPerProc[i], 1);
      j += send_lens_[i];
    }
  } else {
    total_to_send = number_of_vertices_ * 3;

    for (; i < processors_; ++i) {
      send_displs_[i] = j;
      idxIntoSendArray[i] = j;
      send_lens_[i] = localVPerProc[i] * 3;
      j += send_lens_[i];
    }
  }

  send_array_.reserve(total_to_send);

  /* compute the send data_ */

  if (number_of_partitions_ == 0) {
    for (i = 0; i < number_of_vertices_; ++i) {
      j = vToProc[i];
      start_offset = idxIntoSendArray[j];
      send_array_[start_offset] = minimum_vertex_index_ + i;
      send_array_[start_offset + 1] = vertex_weights_[i];
      idxIntoSendArray[j] += 2;
    }
  } else {
    for (i = 0; i < number_of_vertices_; ++i) {
      j = vToProc[i];
      start_offset = idxIntoSendArray[j];
      send_array_[start_offset] = minimum_vertex_index_ + i;
      send_array_[start_offset + 1] = vertex_weights_[i];
      send_array_[start_offset + 2] = partition_vector_[i];
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
  total_to_receive = j;

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
    while (i < total_to_receive) {
      to_origin_vertex_[j] = receive_array_[i++];
      vertex_weights_[j++] = receive_array_[i];
      vertex_weight_ += receive_array_[i++];
    }
  } else {
    while (i < total_to_receive) {
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

  int max_local_vertex = minimum_vertex_index_ + number_of_vertices_;
  int vertices_per_proc_fine = total_number_of_vertices_ / processors_;
  int vToOrigVexist = to_origin_vertex_.capacity();
  int total_to_send;
  int total_to_receive;
  int vertex;
  int start_offset;
  int arrayLen;
  int *array;

  dynamic_array<int> oldIndexToNew(number_of_vertices_);
  dynamic_array<int> minIndexOfMyVerts(processors_);
  dynamic_array<int> minIndexOnProc(processors_);
  dynamic_array<int> idxIntoSendArray(processors_);
  dynamic_array<int> copy_of_requests;

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

  ds::bit_field sent_requests(total_number_of_vertices_);
  sent_requests.unset();

  /* compute all the requests for new indices of remote vertices */

  for (i = 0; i < number_of_pins_; ++i) {
    vertex = pin_list_[i];

    if (vertex < minimum_vertex_index_ || vertex >= max_local_vertex) {
      if (!sent_requests(vertex)) {
        j = std::min(vertex / vertices_per_proc_fine, processors_ - 1);
        data_out_sets_[j]->assign(send_lens_[j]++, vertex);
        sent_requests.set(vertex);
      }
    }
  }

  /* compute number of elements to send to other processors */

  j = 0;
  for (i = 0; i < processors_; ++i) {
    send_displs_[i] = j;
    j += send_lens_[i];
  }

  send_array_.reserve(j);
  copy_of_requests.reserve(j);
  total_to_send = j;

  j = 0;
  for (ij = 0; ij < processors_; ++ij) {
    array = data_out_sets_[ij]->data();
    arrayLen = send_lens_[ij];

    for (i = 0; i < arrayLen; ++i) {
      vertex = array[i];
      send_array_[j] = vertex;
      copy_of_requests[j++] = vertex;
    }
  }

#ifdef DEBUG_HYPERGRAPH
  assert(j == total_to_send);
#endif

  MPI_Alltoall(send_lens_.data(), 1, MPI_INT, receive_lens_.data(), 1, MPI_INT,
               comm);

  /* compute number of elements to receive from other processors */

  j = 0;
  for (i = 0; i < processors_; ++i) {
    receive_displs_[i] = j;
    j += receive_lens_[i];
  }

  receive_array_.reserve(j);
  total_to_receive = j;

  MPI_Alltoallv(send_array_.data(), send_lens_.data(),
                send_displs_.data(), MPI_INT, receive_array_.data(),
                receive_lens_.data(), receive_displs_.data(), MPI_INT, comm);

  /*
    now have received all requests and sent out our requests
    the reply communication will have the dual dimensions
  */

  send_array_.reserve(total_to_receive);

  for (i = 0; i < total_to_receive; ++i) {
#ifdef DEBUG_HYPERGRAPH
    assert(receive_array_[i] >= minVertexIndex &&
           receive_array_[i] < max_local_vertex);
#endif
    send_array_[i] = oldIndexToNew[receive_array_[i] - minimum_vertex_index_];
#ifdef DEBUG_HYPERGRAPH
    assert(send_array_[i] >= 0 && send_array_[i] < numTotalVertices);
#endif
  }

  receive_array_.reserve(total_to_send);

  MPI_Alltoallv(send_array_.data(), receive_lens_.data(),
                receive_displs_.data(), MPI_INT, receive_array_.data(),
                send_lens_.data(), send_displs_.data(), MPI_INT, comm);

  /*
     now the requested vertices are in the copy_of_requests data_
     while their corresponding matchVector values are in
     the corresponding location in the receive_array_
  */

  if (number_of_pins_ < total_number_of_vertices_ / 2) {
    ds::map_from_pos_int<int> stored_requests(number_of_pins_);

    for (i = 0; i < total_to_send; ++i) {
#ifdef DEBUG_HYPERGRAPH
      assert(receive_array_[i] >= 0 && receive_array_[i] < numTotalVertices);
#endif
      stored_requests.insert(copy_of_requests[i], receive_array_[i]);
    }

#ifdef DEBUG_HYPERGRAPH
    for (i = 0; i < numLocalVertices; ++i)
      assert(oldIndexToNew[i] >= 0 && oldIndexToNew[i] < numTotalVertices);
#endif

    /*  now we will convert the local pin list and hash keys */

    for (i = 0; i < number_of_pins_; ++i) {
      vertex = pin_list_[i];

      if (vertex >= minimum_vertex_index_ && vertex < max_local_vertex) {
        pin_list_[i] = oldIndexToNew[vertex - minimum_vertex_index_];
      } else {
        pin_list_[i] = stored_requests.get(vertex);
#ifdef DEBUG_HYPERGRAPH
        assert(localPins[i] >= 0 && localPins[i] < numTotalVertices);
#endif
      }
    }
  } else {
    dynamic_array<int> stored_requests(total_number_of_vertices_);

    for (i = 0; i < total_to_send; ++i) {
#ifdef DEBUG_HYPERGRAPH
      assert(receive_array_[i] >= 0 && receive_array_[i] < numTotalVertices);
#endif
      stored_requests[copy_of_requests[i]] = receive_array_[i];
    }

#ifdef DEBUG_HYPERGRAPH
    for (i = 0; i < numLocalVertices; ++i)
      assert(oldIndexToNew[i] >= 0 && oldIndexToNew[i] < numTotalVertices);
#endif

    /*  now we will convert the local pin list and hash keys */

    for (i = 0; i < number_of_pins_; ++i) {
      vertex = pin_list_[i];

      if (vertex >= minimum_vertex_index_ && vertex < max_local_vertex) {
        pin_list_[i] = oldIndexToNew[vertex - minimum_vertex_index_];
      } else {
        pin_list_[i] = stored_requests[vertex];
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
    total_to_send = number_of_vertices_ * ij;
#ifdef DEBUG_HYPERGRAPH
    assert(total_to_send >= 0);
#endif
    for (i = 0; i < processors_; ++i) {
      send_displs_[i] = j;
      idxIntoSendArray[i] = j;
      send_lens_[i] = localVPerProc[i] * ij;
      j += send_lens_[i];
    }
  } else {
    ij = number_of_partitions_ + 3;
    total_to_send = number_of_vertices_ * ij;
#ifdef DEBUG_HYPERGRAPH
    assert(total_to_send >= 0);
#endif
    for (i = 0; i < processors_; ++i) {
      send_displs_[i] = j;
      idxIntoSendArray[i] = j;
      send_lens_[i] = localVPerProc[i] * ij;
      j += send_lens_[i];
    }
  }

#ifdef DEBUG_HYPERGRAPH
  assert(j == total_to_send);
  assert(total_to_send >= 0);
#endif

  send_array_.reserve(total_to_send);

  /* compute the send data_ */

  if (!vToOrigVexist) {
    for (i = 0; i < number_of_vertices_; ++i) {
      j = vToProc[i];
      start_offset = idxIntoSendArray[j];
      send_array_[start_offset++] = vertex_weights_[i];
      send_array_[start_offset++] = mapToOrigV[i];

      for (ij = 0; ij < number_of_partitions_; ++ij)
        send_array_[start_offset++] =
            partition_vector_[partition_vector_offsets_[ij] + i];

      idxIntoSendArray[j] += (ij + 2);
    }
  } else {
    for (i = 0; i < number_of_vertices_; ++i) {
      j = vToProc[i];
      start_offset = idxIntoSendArray[j];
      send_array_[start_offset++] = to_origin_vertex_[i];
      send_array_[start_offset++] = vertex_weights_[i];
      send_array_[start_offset++] = mapToOrigV[i];

      for (ij = 0; ij < number_of_partitions_; ++ij)
        send_array_[start_offset++] =
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
  total_to_receive = j;

  MPI_Alltoallv(send_array_.data(), send_lens_.data(),
                send_displs_.data(), MPI_INT, receive_array_.data(),
                receive_lens_.data(), receive_displs_.data(), MPI_INT, comm);

  /* change the local vertex information */

  i = 0;
  j = 0;
  vertex_weight_ = 0;

  if (!vToOrigVexist) {
    while (i < total_to_receive) {
      vertex_weights_[j] = receive_array_[i];
      vertex_weight_ += receive_array_[i++];
      mapToOrigV[j] = receive_array_[i++];

      for (ij = 0; ij < number_of_partitions_; ++ij)
        partition_vector_[partition_vector_offsets_[ij] + j] = receive_array_[i++];

      ++j;
    }
  } else {
    to_origin_vertex_.reserve(number_of_vertices_);

    while (i < total_to_receive) {
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

  int max_local_vertex = minimum_vertex_index_ + number_of_vertices_;
  int vertices_per_proc_fine = total_number_of_vertices_ / processors_;
  int vToOrigVexist = to_origin_vertex_.capacity();
  int total_to_send;
  int total_to_receive;
  int vertex;
  int start_offset;
  int arrayLen;
  int *array;

  dynamic_array<int> oldIndexToNew(number_of_vertices_);
  dynamic_array<int> minIndexOfMyVerts(processors_);
  dynamic_array<int> minIndexOnProc(processors_);
  dynamic_array<int> idxIntoSendArray(processors_);
  dynamic_array<int> copy_of_requests;

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

  ds::bit_field sent_requests(total_number_of_vertices_);
  sent_requests.unset();

  /* compute all the requests for new indices of remote vertices  */

  for (i = 0; i < number_of_pins_; ++i) {
    vertex = pin_list_[i];

    if (vertex < minimum_vertex_index_ || vertex >= max_local_vertex) {
      if (!sent_requests(vertex)) {
        j = std::min(vertex / vertices_per_proc_fine, processors_ - 1);
        data_out_sets_[j]->assign(send_lens_[j]++, vertex);
        sent_requests.set(vertex);
      }
    }
  }

  /* compute number of elements to send to other processors */

  j = 0;
  for (i = 0; i < processors_; ++i) {
    send_displs_[i] = j;
    j += send_lens_[i];
  }

  send_array_.reserve(j);
  copy_of_requests.reserve(j);
  total_to_send = j;

  j = 0;
  for (ij = 0; ij < processors_; ++ij) {
    array = data_out_sets_[ij]->data();
    arrayLen = send_lens_[ij];

    for (i = 0; i < arrayLen; ++i) {
      vertex = array[i];
      send_array_[j] = vertex;
      copy_of_requests[j++] = vertex;
    }
  }

#ifdef DEBUG_HYPERGRAPH
  assert(j == total_to_send);
#endif

  MPI_Alltoall(send_lens_.data(), 1, MPI_INT, receive_lens_.data(), 1, MPI_INT,
               comm);

  /* compute number of elements to receive from other processors */

  j = 0;
  for (i = 0; i < processors_; ++i) {
    receive_displs_[i] = j;
    j += receive_lens_[i];
  }

  receive_array_.reserve(j);
  total_to_receive = j;

  MPI_Alltoallv(send_array_.data(), send_lens_.data(),
                send_displs_.data(), MPI_INT, receive_array_.data(),
                receive_lens_.data(), receive_displs_.data(), MPI_INT, comm);

  /*
    now have received all requests and sent out our requests
    the reply communication will have the dual dimensions
  */

  send_array_.reserve(total_to_receive);

  for (i = 0; i < total_to_receive; ++i) {
#ifdef DEBUG_HYPERGRAPH
    assert(receive_array_[i] >= minVertexIndex &&
           receive_array_[i] < max_local_vertex);
#endif
    send_array_[i] = oldIndexToNew[receive_array_[i] - minimum_vertex_index_];
#ifdef DEBUG_HYPERGRAPH
    assert(send_array_[i] >= 0 && send_array_[i] < numTotalVertices);
#endif
  }

  receive_array_.reserve(total_to_send);

  MPI_Alltoallv(send_array_.data(), receive_lens_.data(),
                receive_displs_.data(), MPI_INT, receive_array_.data(),
                send_lens_.data(), send_displs_.data(), MPI_INT, comm);

  /*
  now the requested vertices are in the copy_of_requests data_
   while their corresponding matchVector values are in
   the corresponding location in the receive_array_
  */

  if (number_of_pins_ < total_number_of_vertices_ / 2) {
    ds::map_from_pos_int<int> stored_requests(number_of_pins_);

    for (i = 0; i < total_to_send; ++i) {
#ifdef DEBUG_HYPERGRAPH
      assert(receive_array_[i] >= 0 && receive_array_[i] < numTotalVertices);
#endif
      stored_requests.insert(copy_of_requests[i], receive_array_[i]);
    }

/* now we will convert the local pin list and hash keys */

#ifdef DEBUG_HYPERGRAPH
    for (i = 0; i < numLocalVertices; ++i)
      assert(oldIndexToNew[i] >= 0 && oldIndexToNew[i] < numTotalVertices);
#endif

    for (i = 0; i < number_of_pins_; ++i) {
      vertex = pin_list_[i];

      if (vertex >= minimum_vertex_index_ && vertex < max_local_vertex) {
        pin_list_[i] = oldIndexToNew[vertex - minimum_vertex_index_];
      } else {
        pin_list_[i] = stored_requests.get(vertex);
#ifdef DEBUG_HYPERGRAPH
        assert(localPins[i] >= 0 && localPins[i] < numTotalVertices);
#endif
      }
    }
  } else {
    dynamic_array<int> stored_requests(total_number_of_vertices_);

    for (i = 0; i < total_to_send; ++i) {
#ifdef DEBUG_HYPERGRAPH
      assert(receive_array_[i] >= 0 && receive_array_[i] < numTotalVertices);
#endif
      stored_requests[copy_of_requests[i]] = receive_array_[i];
    }

/* now we will convert the local pin list and hash keys */

#ifdef DEBUG_HYPERGRAPH
    for (i = 0; i < numLocalVertices; ++i)
      assert(oldIndexToNew[i] >= 0 && oldIndexToNew[i] < numTotalVertices);
#endif

    for (i = 0; i < number_of_pins_; ++i) {
      vertex = pin_list_[i];

      if (vertex >= minimum_vertex_index_ && vertex < max_local_vertex) {
        pin_list_[i] = oldIndexToNew[vertex - minimum_vertex_index_];
      } else {
        pin_list_[i] = stored_requests[vertex];
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
    total_to_send = number_of_vertices_ * ij;

    for (i = 0; i < processors_; ++i) {
      send_displs_[i] = j;
      idxIntoSendArray[i] = j;
      send_lens_[i] = localVPerProc[i] * ij;
      j += send_lens_[i];
    }
  } else {
    ij = number_of_partitions_ + 4;
    total_to_send = number_of_vertices_ * ij;

    for (i = 0; i < processors_; ++i) {
      send_displs_[i] = j;
      idxIntoSendArray[i] = j;
      send_lens_[i] = localVPerProc[i] * ij;
      j += send_lens_[i];
    }
  }

  send_array_.reserve(total_to_send);

  /* compute the send data_ */

  if (!vToOrigVexist) {
    for (i = 0; i < number_of_vertices_; ++i) {
      j = vToProc[i];
      start_offset = idxIntoSendArray[j];
      send_array_[start_offset++] = vertex_weights_[i];
      send_array_[start_offset++] = mapToInterV[i];
      send_array_[start_offset++] = mapToOrigV[i];

      for (ij = 0; ij < number_of_partitions_; ++ij)
        send_array_[start_offset++] =
            partition_vector_[partition_vector_offsets_[ij] + i];

      idxIntoSendArray[j] += (3 + ij);
    }
  } else {
    for (i = 0; i < number_of_vertices_; ++i) {
      j = vToProc[i];
      start_offset = idxIntoSendArray[j];
      send_array_[start_offset++] = to_origin_vertex_[i];
      send_array_[start_offset++] = vertex_weights_[i];
      send_array_[start_offset++] = mapToInterV[i];
      send_array_[start_offset++] = mapToOrigV[i];

      for (ij = 0; ij < number_of_partitions_; ++ij)
        send_array_[start_offset++] =
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
  total_to_receive = j;

  MPI_Alltoallv(send_array_.data(), send_lens_.data(),
                send_displs_.data(), MPI_INT, receive_array_.data(),
                receive_lens_.data(), receive_displs_.data(), MPI_INT, comm);

  /* change the local vertex information */

  i = 0;
  j = 0;
  vertex_weight_ = 0;

  if (!vToOrigVexist) {
    while (i < total_to_receive) {
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

    while (i < total_to_receive) {
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

  int max_local_vertex = minimum_vertex_index_ + number_of_vertices_;
  int vertices_per_proc_fine = total_number_of_vertices_ / processors_;
  int total_to_send;
  int total_to_receive;
  int vertex;
  int start_offset;
  int arrayLen;
  int vToOrigVexist = to_origin_vertex_.capacity();
  int *array;

  dynamic_array<int> oldIndexToNew(number_of_vertices_);
  dynamic_array<int> minIndexOfMyVerts(processors_);
  dynamic_array<int> minIndexOnProc(processors_);
  dynamic_array<int> idxIntoSendArray(processors_);
  dynamic_array<int> copy_of_requests;

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

  ds::bit_field sent_requests(total_number_of_vertices_);
  sent_requests.unset();

  /*
     compute all the requests for new indices of remote vertices
     because we need to modify the match vector of the finer graph,
     this needs to include those vertices that are in this match
     vector as well
  */

  for (i = 0; i < number_of_pins_; ++i) {
    vertex = pin_list_[i];

    if (vertex < minimum_vertex_index_ || vertex >= max_local_vertex) {
      if (!sent_requests(vertex)) {
        j = std::min(vertex / vertices_per_proc_fine, processors_ - 1);
        data_out_sets_[j]->assign(send_lens_[j]++, vertex);
        sent_requests.set(vertex);
      }
    }
  }

  for (i = 0; i < numLocalFineVertices; ++i) {
    vertex = fineMatchVector[i];

    if (vertex < minimum_vertex_index_ || vertex >= max_local_vertex) {
      if (!sent_requests(vertex)) {
        j = std::min(vertex / vertices_per_proc_fine, processors_ - 1);
        data_out_sets_[j]->assign(send_lens_[j]++, vertex);
        sent_requests.set(vertex);
      }
    }
  }

  /* compute number of elements to send to other processors */

  j = 0;
  for (i = 0; i < processors_; ++i) {
    send_displs_[i] = j;
    j += send_lens_[i];
  }

  send_array_.reserve(j);
  copy_of_requests.reserve(j);
  total_to_send = j;

  j = 0;
  for (ij = 0; ij < processors_; ++ij) {
    array = data_out_sets_[ij]->data();
    arrayLen = send_lens_[ij];

    for (i = 0; i < arrayLen; ++i) {
      vertex = array[i];
      send_array_[j] = vertex;
      copy_of_requests[j++] = vertex;
    }
  }

#ifdef DEBUG_HYPERGRAPH
  assert(j == total_to_send);
#endif

  MPI_Alltoall(send_lens_.data(), 1, MPI_INT, receive_lens_.data(), 1, MPI_INT,
               comm);

  /* compute number of elements to receive from other processors */

  j = 0;
  for (i = 0; i < processors_; ++i) {
    receive_displs_[i] = j;
    j += receive_lens_[i];
  }

  receive_array_.reserve(j);
  total_to_receive = j;

  MPI_Alltoallv(send_array_.data(), send_lens_.data(),
                send_displs_.data(), MPI_INT, receive_array_.data(),
                receive_lens_.data(), receive_displs_.data(), MPI_INT, comm);

  /*
    now have received all requests and sent out our requests
    the reply communication will have the dual dimensions
  */

  send_array_.reserve(total_to_receive);

  for (i = 0; i < total_to_receive; ++i) {
#ifdef DEBUG_HYPERGRAPH
    assert(receive_array_[i] >= minVertexIndex &&
           receive_array_[i] < max_local_vertex);
#endif
    send_array_[i] = oldIndexToNew[receive_array_[i] - minimum_vertex_index_];
#ifdef DEBUG_HYPERGRAPH
    assert(send_array_[i] >= 0 && send_array_[i] < numTotalVertices);
#endif
  }

  receive_array_.reserve(total_to_send);

  MPI_Alltoallv(send_array_.data(), receive_lens_.data(),
                receive_displs_.data(), MPI_INT, receive_array_.data(),
                send_lens_.data(), send_displs_.data(), MPI_INT, comm);

  /*
    now the requested vertices are in the copy_of_requests data_
    while their corresponding matchVector values are in
  */

  if (number_of_pins_ < total_number_of_vertices_ / 2) {
    ds::map_from_pos_int<int> stored_requests(number_of_pins_);

    for (i = 0; i < total_to_send; ++i) {
#ifdef DEBUG_HYPERGRAPH
      assert(receive_array_[i] >= 0 && receive_array_[i] < numTotalVertices);
#endif
      stored_requests.insert(copy_of_requests[i], receive_array_[i]);
    }

/* now we will convert the local pin list and hash keys */

#ifdef DEBUG_HYPERGRAPH
    for (i = 0; i < numLocalVertices; ++i)
      assert(oldIndexToNew[i] >= 0 && oldIndexToNew[i] < numTotalVertices);
#endif

    for (i = 0; i < number_of_pins_; ++i) {
      vertex = pin_list_[i];

      if (vertex >= minimum_vertex_index_ && vertex < max_local_vertex) {
        pin_list_[i] = oldIndexToNew[vertex - minimum_vertex_index_];
      } else {
        pin_list_[i] = stored_requests.get(vertex);
#ifdef DEBUG_HYPERGRAPH
        assert(localPins[i] >= 0 && localPins[i] < numTotalVertices);
#endif
      }
    }

    /* now we will convert the fine graph match vector */

    for (i = 0; i < numLocalFineVertices; ++i) {
      vertex = fineMatchVector[i];

      if (vertex >= minimum_vertex_index_ && vertex < max_local_vertex) {
        fineMatchVector[i] = oldIndexToNew[vertex - minimum_vertex_index_];
      } else {
        fineMatchVector[i] = stored_requests.get(vertex);
      }
#ifdef DEBUG_HYPERGRAPH
      assert(fineMatchVector[i] >= 0 && fineMatchVector[i] < numTotalVertices);
#endif
    }
  } else {
    dynamic_array<int> stored_requests(total_number_of_vertices_);

    for (i = 0; i < total_to_send; ++i) {
#ifdef DEBUG_HYPERGRAPH
      assert(receive_array_[i] >= 0 && receive_array_[i] < numTotalVertices);
#endif
      stored_requests[copy_of_requests[i]] = receive_array_[i];
    }

/* now we will convert the local pin list and hash keys */

#ifdef DEBUG_HYPERGRAPH
    for (i = 0; i < numLocalVertices; ++i)
      assert(oldIndexToNew[i] >= 0 && oldIndexToNew[i] < numTotalVertices);
#endif

    for (i = 0; i < number_of_pins_; ++i) {
      vertex = pin_list_[i];

      if (vertex >= minimum_vertex_index_ && vertex < max_local_vertex) {
        pin_list_[i] = oldIndexToNew[vertex - minimum_vertex_index_];
      } else {
        pin_list_[i] = stored_requests[vertex];
#ifdef DEBUG_HYPERGRAPH
        assert(localPins[i] >= 0 && localPins[i] < numTotalVertices);
#endif
      }
    }

    /* now we will convert the fine graph match vector */

    for (i = 0; i < numLocalFineVertices; ++i) {
      vertex = fineMatchVector[i];

      if (vertex >= minimum_vertex_index_ && vertex < max_local_vertex) {
        fineMatchVector[i] = oldIndexToNew[vertex - minimum_vertex_index_];
      } else {
        fineMatchVector[i] = stored_requests[vertex];
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
    total_to_send = number_of_vertices_ * ij;

    for (i = 0; i < processors_; ++i) {
      send_displs_[i] = j;
      idxIntoSendArray[i] = j;
      send_lens_[i] = localVPerProc[i] * ij;
      j += send_lens_[i];
    }
  } else {
    ij = number_of_partitions_ + 1;
    total_to_send = number_of_vertices_ * ij;

    for (i = 0; i < processors_; ++i) {
      send_displs_[i] = j;
      idxIntoSendArray[i] = j;
      send_lens_[i] = localVPerProc[i] * ij;
      j += send_lens_[i];
    }
  }

#ifdef DEBUG_HYPERGRAPH
  assert(j == total_to_send);
#endif

  send_array_.reserve(total_to_send);

  /* compute the send data_ */

  if (vToOrigVexist > 0) {
    for (i = 0; i < number_of_vertices_; ++i) {
      j = vToProc[i];

      start_offset = idxIntoSendArray[j];
      send_array_[start_offset++] = to_origin_vertex_[i];
      send_array_[start_offset++] = vertex_weights_[i];

      for (ij = 0; ij < number_of_partitions_; ++ij)
        send_array_[start_offset++] =
            partition_vector_[partition_vector_offsets_[ij] + i];

      idxIntoSendArray[j] += (ij + 2);
    }
  } else {
    for (i = 0; i < number_of_vertices_; ++i) {
      j = vToProc[i];

      start_offset = idxIntoSendArray[j];
      send_array_[start_offset++] = vertex_weights_[i];

      for (ij = 0; ij < number_of_partitions_; ++ij)
        send_array_[start_offset++] =
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
  total_to_receive = j;

  MPI_Alltoallv(send_array_.data(), send_lens_.data(),
                send_displs_.data(), MPI_INT, receive_array_.data(),
                receive_lens_.data(), receive_displs_.data(), MPI_INT, comm);

  /* change the local vertex information */

  i = 0;
  j = 0;
  vertex_weight_ = 0;

  if (vToOrigVexist > 0) {
    while (i < total_to_receive) {
      to_origin_vertex_[j] = receive_array_[i++];
      vertex_weights_[j] = receive_array_[i];
      vertex_weight_ += receive_array_[i++];

      for (ij = 0; ij < number_of_partitions_; ++ij)
        partition_vector_[partition_vector_offsets_[ij] + j] = receive_array_[i++];

      ++j;
    }
  } else {
    while (i < total_to_receive) {
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
  int max_local_vertex = minimum_vertex_index_ + number_of_vertices_;
  int vPerProc = total_number_of_vertices_ / processors_;
  int locCutsize = 0;
  int totCutsize;
  int numSpanned;
  int totToSend;
  int totToRecv;
  int vertex;
  int end_offset;
  int part;

  int *pVector = &partition_vector_[partition_vector_offsets_[pNum]];
  int *array;

  dynamic_array<int> spannedPart(numParts);
  dynamic_array<int> copy_of_requests;

  send_lens_.assign(processors_, 0);

  ds::bit_field sent_requests(total_number_of_vertices_);
  sent_requests.unset();

  // Compute all the requests for partition vector values of remote vertices.
  assert(number_of_pins_ == pin_list_.size());
  for (const auto &pin : pin_list_) {
    if ((pin < minimum_vertex_index_ || pin >= max_local_vertex)
        && !sent_requests[pin]) {
      int j = std::min(pin / vPerProc, processors_ - 1);
      data_out_sets_[j][send_lens_[j]++] = pin;
      sent_requests.set(pin);
    }
  }

  // Compute number of elements to send to other processors.
  int total_to_send = 0;
  for (i = 0; i < processors_; ++i) {
    send_displs_[i] = total_to_send;
    total_to_send += send_lens_[i];
  }

  send_array_.assign(total_to_send, 0);
  copy_of_requests.assign(total_to_send, 0;

  int ij = 0;
  for (int i = 0; i < processors_; ++i) {
    auto array = data_out_sets_[i];
    int array_len = send_lens_[i];

    for (int j = 0; j < array_len; ++j) {
      vertex = array[j];
      send_array_[ij] = vertex;
      copy_of_requests[ij++] = vertex;
    }
  }

  MPI_Alltoall(send_lens_.data(), 1, MPI_INT, receive_lens_.data(), 1, MPI_INT,
               comm);

  // Compute number of elements to receive from other processors
  ij = 0;
  for (int i = 0; i < processors_; ++i) {
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
           receive_array_[i] < max_local_vertex);
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
  // now the requested vertices are in the copy_of_requests data_
  // while their corresponding matchVector values are in
  // the corresponding location in the receive_array_
  // ###

  if (number_of_pins_ < total_number_of_vertices_ / 2) {
    ds::map_from_pos_int<int> stored_requests(number_of_pins_);

    for (i = 0; i < totToSend; ++i) {
#ifdef DEBUG_HYPERGRAPH
      assert(receive_array_[i] >= 0 && receive_array_[i] < numParts);
#endif
      stored_requests.insert(copy_of_requests[i], receive_array_[i]);
    }

    // ###
    // now procees to compute the contribution to
    // cutsize of the locally held hyperedges
    // ###

    for (i = 0; i < number_of_hyperedges_; ++i) {
      for (j = 0; j < numParts; ++j)
        spannedPart[j] = 0;

      numSpanned = 0;
      end_offset = hyperedge_offsets_[i + 1];

      for (j = hyperedge_offsets_[i]; j < end_offset; ++j) {
        ij = pin_list_[j];

        if (ij < minimum_vertex_index_ || ij >= max_local_vertex) {
          part = stored_requests.get(ij);
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
    dynamic_array<int> stored_requests(total_number_of_vertices_);

    for (i = 0; i < totToSend; ++i) {
#ifdef DEBUG_HYPERGRAPH
      assert(receive_array_[i] >= 0 && receive_array_[i] < numParts);
#endif
      stored_requests[copy_of_requests[i]] = receive_array_[i];
    }

    // ###
    // now procees to compute the contribution to
    // cutsize of the locally held hyperedges
    // ###

    for (i = 0; i < number_of_hyperedges_; ++i) {
      for (j = 0; j < numParts; ++j)
        spannedPart[j] = 0;

      numSpanned = 0;
      end_offset = hyperedge_offsets_[i + 1];

      for (j = hyperedge_offsets_[i]; j < end_offset; ++j) {
        ij = pin_list_[j];

        if (ij < minimum_vertex_index_ || ij >= max_local_vertex) {
          part = stored_requests[ij];
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


void hypergraph::check_vertex_and_hyperedge_lengths(
    int hyperedge_data_length, const dynamic_array<int> &hyperedge_data,
    MPI_Comm comm) {
  int hyperedges_in_file = 0;
  int j = 0;
  for (int i = 0; i < hyperedge_data_length; i += hyperedge_data[i]) {
    int hyperedge_index = hyperedge_data[i] - 2;
    if (hyperedge_index > j) {
      j = hyperedge_index;
    }
    ++hyperedges_in_file;
  }

  int max_hyperedge_length;
  int number_of_edges;
  MPI_Allreduce(&j, &max_hyperedge_length, 1, MPI_INT, MPI_MAX, comm);
  MPI_Reduce(&hyperedges_in_file, &number_of_edges, 1, MPI_INT, MPI_SUM, 0,
             comm);

  Funct::setMaxHedgeLen(max_hyperedge_length);

  if (rank_ == 0 && display_option_ > 0) {
    out << "|--- Hypergraph " << filename << " (on file):" << std::endl
        << "| |V| = " << total_number_of_vertices_  << std::endl;
        << "| |E| = " << number_of_edges << std::endl;
  }
}


void hypergraph::load_data_from_blocks(
    const int data_length, const dynamic_array<int> &hypergraph_data) {
  // Hyperedges stored on file in blocks:
  // [0]  =  block capacity
  // [1]  =  hyperedge weight
  // [2]-[block capacity-1] = hyperedge vertices
  //
  int i = 0;
  int hyperedge_index = 0;
  int pin_counter = 0;
  while (i < hyperedge_data_length) {
    int chunk = hyperedge_data[i];
    int length = chunk - 2;

    if (length > 1) {
      hyperedge_weights_[hyperedge_index] = hyperedge_data[i + 1];
      hyperedge_offsets_[hyperedge_index++] = pin_counter;
      for (int j = 2; j < chunk; ++j) {
        pin_list_[pin_counter++] = hyperedge_data[i + j];
      }
    }
    i += chunk;
  }

  hyperedge_offsets_[hyperedge_index] = pin_counter;

  number_of_pins_ = pin_counter;
  number_of_hyperedges_ = hyperedge_index;
  do_not_coarsen = 0;
  number_of_partitions_ = 0;
}


void hypergraph::check_loaded_vertex_and_hyperedge_lengths() const {
  if (display_option_ > 0) {
    int pins;
    int edges;
    MPI_Reduce(&number_of_pins_, &pins, 1, MPI_INT, MPI_SUM, 0, comm);
    MPI_Reduce(&number_of_hyperedges_, &edges, 1, MPI_INT, MPI_SUM, 0, comm);
    if (rank_ == 0) {
      out << "|--- Hypergraph " << filename << " (as loaded):" << std::endl
          << "| |V| = " << total_number_of_vertices_  << std::endl
          << "| |E| = " << edges << std::endl
          << "| |Pins| = " << pins << std::endl
          << "| # Processors = " << processors_ << std::endl
          << "| " << std::endl;
    }
  }
}

int hypergraph::compute_number_of_elements_to_send(
    dynamic_array<int> &copy_of_requests) {
  int total_to_send = 0;
  for (int i = 0; i < processors_; ++i) {
    send_displs_[i] = total_to_send;
    total_to_send += send_lens_[i];
  }

  send_array_.resize(total_to_send);
  copy_of_requests.resize(total_to_send);

  int j = 0;
  for (int p = 0; p < processors_; ++p) {
    for (int i = 0; i < send_lens_[p]; ++i) {
      vertex = data_out_sets_[p][i];
      send_array_[j] = vertex;
      copy_of_requests[j++] = vertex;
    }
  }

  return total_to_send;
}

void hypergraph::compute_number_of_elements_to_receive() {
  int to_receive = 0;
  for (i = 0; i < processors_; ++i) {
    receive_displs_[i] = to_receive;
    to_receive += receive_lens_[i];
  }

  receive_array_.resize(total_to_receive);
  return to_receive;
}

int hypergraph::get_processor(int vertex, int vertices_per_processor,
                              ds::complete_binary_tree<int> *vertex_to_proc) {
  if (vertex_to_proc == nullptr) {
    return std::min(vertex / vertices_per_processor, processors_ - 1);
  } else {
    return vertex_to_proc->root_value(vertex);
  }
}


void hypergraph::compute_requests_for_remote_vertex_matches(
    int vertices_per_processor, ds::complete_binary_tree<int> *vertex_to_proc) {
  send_lens_.assign(processors_, 0);
  std::bitset<total_number_of_vertices_> sent_requests;

  // Compute all the requests for remote vertex matches
  for (int i = 0; i < number_of_pins_; ++i) {
    int vertex = pin_list_[i];
    if (vertex < minimum_vertex_index_ || vertex >= max_local_vertex) {
      if (!sent_requests.test(vertex)) {
        int p = get_processor(vertex, vertices_per_processor, vertex_to_proc);
        data_out_sets_[p][send_lens_[p]++] = vertex;
        sent_requests.set(vertex);
      }
      original_contracted_pin_list[i] = -1;
    } else {
      original_contracted_pin_list[i] =
          match_vector_[vertex - minimum_vertex_index_];
    }
  }
}

void hypergraph::choose_non_local_vertices_format(
    dynamic_array<int> &original_contracted_pin_list) {
  if (number_of_pins_ < total_number_of_vertices_ / 2) {
    ds::map_from_pos_int<int> stored_requests(number_of_pins_);
    for (int i = 0; i < total_to_send; ++i) {
      stored_requests.insert(copy_of_requests[i], receive_array_[i]);
    }
    // contract remaining local pins
    for (int i = 0; i < number_of_pins_; ++i) {
      if (original_contracted_pin_list[i] == -1) {
        original_contracted_pin_list[i] = stored_requests.get(pin_list_[i]);
      }
    }
  } else {
    dynamic_array<int> non_local_matches(total_number_of_vertices_);
    for (int i = 0; i < total_to_send; ++i) {
      non_local_matches[copy_of_requests[i]] = receive_array_[i];
    }
    // contract remaining local pins
    for (int i = 0; i < number_of_pins_; ++i) {
      if (original_contracted_pin_list[i] == -1) {
        original_contracted_pin_list[i] = non_local_matches[pin_list_[i]];
      }
    }
  }
}

void hypergraph::send_coarse_hyperedges(
    ds::dynamic_array<int> &original_contracted_pin_list,
    int &total_to_send, int &total_to_receive) {
  for (int i = 0; i < number_of_hyperedges_; ++i) {
    int j = hyperedge_offsets_[i + 1] - hyperedge_offsets_[i];
    int start = hyperedge_offsets_[i];
    original_contracted_pin_list.sort_between(start, start + j - 1);
  }

  int contracted_pin_list_length = 0;
  int max_local_coarse_hedge_length = 0;
  int number_of_contracted_hedges = 0;
  int contracted_hedge_length = 0;
  dynamic_array<int> contracted_hedge_offsets(number_of_hyperedges_);
  dynamic_array<int> contracted_hedge_weights(number_of_hyperedges_);
  dynamic_array<int> contracted_pin_list;
  for (int i = 0; i < number_of_hyperedges_; ++i) {
    contracted_hedge_offsets[number_of_contracted_hedges] = contracted_pin_list_length;
    int end_offset = hyperedge_offsets_[i + 1];
    int start_offset = hyperedge_offsets_[i];

    contracted_pin_list[contracted_pin_list_length] =
        original_contracted_pin_list[start_offset];
    contracted_hedge_length = 1;

    for (int j = start_offset + 1; j < end_offset; ++j) {
      int index = contracted_pin_list_length + contracted_hedge_length;
      if (original_contracted_pin_list[j] != contracted_pin_list[index - 1]) {
        contracted_pin_list[index] = original_contracted_pin_list[j];
        contracted_hedge_length++;
      }
    }

    if (contracted_hedge_length > 1) {
      contracted_pin_list_length += contracted_hedge_length;
      contracted_hedge_weights[number_of_contracted_hedges++] = hyperedge_weights_[i];
      if (contracted_hedge_length > max_local_coarse_hedge_length) {
        max_local_coarse_hedge_length = contracted_hedge_length;
      }
    }
  }

  contracted_hedge_offsets[number_of_contracted_hedges] = contracted_pin_list_length;

  // - compute max_coarse_hedge_length and set in Funct
  // - compute hash-keys for each hyperedge
  // - send hyperedges to processor determined by corresponding hash key
  int max_coarse_hedge_length;
  MPI_Allreduce(&max_local_coarse_hedge_length, &max_coarse_hedge_length, 1,
                MPI_INT, MPI_MAX, comm);

  Funct::setMaxHedgeLen(max_coarse_hedge_length);

  send_lens_.assign(processors_, 0);
  for (int i = 0; i < number_of_contracted_hedges; ++i) {
    int start_offset = contracted_hedge_offsets[i];
    int end_offset = contracted_hedge_offsets[i + 1];
    int contracted_hedge_length = end_offset - start_offset;

    int p = Funct::computeHash(&contracted_pin_list[start_offset],
                               contracted_hedge_length) % processors_;

    data_out_sets_[p][send_lens_[p]++] = contracted_hedge_length + 2;
    data_out_sets_[p][send_lens_[p]++] = contracted_hedge_weights[i];

    for (int j = start_offset; j < end_offset; ++j) {
      data_out_sets_[p][send_lens_[p]++] = contracted_pin_list[j];
    }
  }

  // Compute number of elements to send to other processors
  total_to_send = compute_number_of_elements_to_send(copy_of_requests);
  MPI_Alltoall(send_lens_.data(), 1, MPI_INT, receive_lens_.data(), 1, MPI_INT,
               comm);

  // Compute number of elements to receive from other processors.
  total_to_receive = compute_number_of_elements_to_receive();
  MPI_Alltoallv(send_array_.data(), send_lens_.data(),
                send_displs_.data(), MPI_INT, receive_array_.data(),
                receive_lens_.data(), receive_displs_.data(), MPI_INT, comm);
}

void hypergraph::process_new_hyperedges(hypergraph &coarse,
                                        ds::new_hyperedge_index_table &table) {
  dynamic_array<int> coarse_local_pins;
  dynamic_array<int> coarse_hedge_offsets;
  dynamic_array<int> coarse_hedge_weights;
  int number_coarse_pins = 0;
  int number_course_hedges = 0;
  coarse_hedge_offsets[number_course_hedges] = number_coarse_pins;

  int i = 0;
  while (i < total_to_receive) {
    int coarse_hyperedge_length = receive_array_[i] - 2;
    HashKey hash_key = Funct::computeHash(&receive_array_[i + 2],
                                          coarse_hyperedge_length);

    // See if a duplicate hyperedge exists
    int number_seen = -1;
    int duplicated_hyperedge = -1;
    do {
      int try_hyperedge = table.getHedgeIndex(hash_key, number_seen);
      if (try_hyperedge >= 0) {
        int start_offset = coarse_hedge_offsets[try_hyperedge];
        int end_offset = coarse_hedge_offsets[try_hyperedge + 1];
        int length = end_offset - start_offset;
        if (length == coarse_hyperedge_length) {
          int j = 0;
          for (j = 0; j < length; ++j) {
            if (receive_array_[i + 2 + j] != coarse_local_pins[start_offset + j]) {
              break;
            }
          }
          if (j == length) {
            duplicated_hyperedge = try_hyperedge;
            break;
          }
        }
      }
    } while (number_seen >= 0);

    if (duplicated_hyperedge == -1) {
      table.insertKey(hash_key, number_course_hedges);
      coarse_hedge_weights[number_course_hedges++] = receive_array_[i + 1];

      for (int j = 0; j < coarse_hyperedge_length; ++j) {
        coarse_local_pins[number_coarse_pins + j] = receive_array_[i + 2 + j];
      }

      number_coarse_pins += coarse_hyperedge_length;
      coarse_hedge_offsets[number_course_hedges] = number_coarse_pins;
    } else {
      coarse_hedge_weights[duplicated_hyperedge] += receive_array_[i + 1];
    }
    i += receive_array_[i];
  }

  // now set the coarse hypergraph
  coarse.set_number_of_hyperedges(number_course_hedges);
  coarse.set_number_of_pins(number_coarse_pins);
  coarse.allocate_hyperedge_memory(number_course_hedges, number_coarse_pins);

  coarse.hyperedge_offsets().set_data(coarse_hedge_offsets.data(), number_course_hedges);
  coarse.hyperedge_weights().set_data(coarse_hedge_weights.data(), number_course_hedges);
  coarse.pin_list().set_data(coarse_local_pins.data(), number_coarse_pins);
}

}  // namespace parallel
}  // namespace parkway
