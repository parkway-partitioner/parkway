// ### ParaGreedyKwayRefiner.cpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 3/2/2005: Last Modified
//
// ###

#include "refiners/parallel/k_way_greedy_refiner.hpp"
#include "utility/random.hpp"
#include "utility/logging.hpp"
#include <iostream>
#include <cmath>

namespace parkway {
namespace parallel {

k_way_greedy_refiner::k_way_greedy_refiner(int rank, int nProcs, int nParts, int
                                           numVperP, int eExit, double lim)
    : refiner(rank, nProcs, nParts) {
  int i;
  int j;
  int ij;

  // ###
  // create the move set structures
  // ###

  early_exit_ = eExit;
  limit_ = lim;
  total_number_of_vertices_moved_ = 0;
  ij = number_of_parts_ * number_of_parts_;

  move_sets_.resize(ij);
  move_set_data_.resize(ij << 1);
  index_into_move_set_.resize(ij);
  number_of_vertices_moved_.resize(ij);

  for (i = 0; i < ij; ++i)
    move_sets_[i] = new dynamic_array<int>;

  for (i = 0; i < number_of_parts_; ++i) {
    ij = i * number_of_parts_;

    for (j = 0; j < number_of_parts_; ++j) {
      if (i == j) {
        index_into_move_set_[ij + j] = -1;
      } else {
        index_into_move_set_[ij + j] = (ij + j) << 1;
      }
    }
  }

  // ###
  // build tables
  // ###

  movement_sets_ = new ds::movement_set_table(number_of_parts_, processors_);

  number_of_neighbor_parts_.resize(0);
  neighbors_of_vertices_.resize(0);
  neighbors_of_vertices_offsets_.resize(0);
  hyperedge_vertices_in_part_.resize(0);
  hyperedge_vertices_in_part_offsets_.resize(0);
  vertices_.resize(0);
  moved_vertices_.resize(0);
  number_of_parts_spanned_.resize(0);
  spanned_parts_.resize(0);
  locked_.reserve(0);
  vertices_seen_.reserve(0);
}

k_way_greedy_refiner::~k_way_greedy_refiner() {
}

void k_way_greedy_refiner::display_options() const {
  info("|--- PARA_REF: \n"
       "|- PKWAY: eeL = %.2f eExit = %i\n|\n", limit_, early_exit_);
}

void k_way_greedy_refiner::release_memory() {
  hyperedge_weights_.resize(0);
  hyperedge_offsets_.resize(0);
  local_pin_list_.resize(0);

  vertex_to_hyperedges_offset_.resize(0);
  vertex_to_hyperedges_.resize(0);
  allocated_hyperedges_.resize(0);

  number_of_neighbor_parts_.resize(0);
  neighbors_of_vertices_.resize(0);
  neighbors_of_vertices_offsets_.resize(0);
  hyperedge_vertices_in_part_.resize(0);
  hyperedge_vertices_in_part_offsets_.resize(0);
  vertices_.resize(0);
  moved_vertices_.resize(0);
  seen_vertices_.resize(0);
  number_of_parts_spanned_.resize(0);
  spanned_parts_.resize(0);

  non_local_vertices_.resize(0);
  part_indices_.resize(0);
  index_into_part_indices_.resize(0);
  non_local_vertices_to_hyperedges_.resize(0);
  non_local_vertices_to_hyperedges_offsets_.resize(0);

  to_non_local_vertices_.destroy();

  free_memory();
}

void k_way_greedy_refiner::initialize_data_structures(
    const parallel::hypergraph &h, MPI_Comm comm) {
  int i;
  int j;
  int ij;

  initialize_partition_structures(h, comm);

  movement_sets_->set_max_part_weight(maximum_part_weight_);

  // ###
  // init data_ structures used in refinement
  // ###

  vertices_seen_.reserve(number_of_local_vertices_);
  vertices_seen_.unset();

  locked_.reserve(number_of_local_vertices_);
  vertices_.resize(number_of_local_vertices_);
  seen_vertices_.resize(number_of_local_vertices_);
  spanned_parts_.resize(number_of_parts_);
  number_of_parts_spanned_.resize(number_of_hyperedges_);
  number_of_neighbor_parts_.resize(number_of_local_vertices_);
  neighbors_of_vertices_.resize(number_of_local_vertices_ * number_of_parts_);
  neighbors_of_vertices_offsets_.resize(number_of_local_vertices_ + 1);

  hyperedge_vertices_in_part_.resize(number_of_hyperedges_ * number_of_parts_);
  hyperedge_vertices_in_part_offsets_.resize(number_of_hyperedges_ + 1);

  j = 0;
  ij = number_of_parts_;

  for (i = 0; i <= number_of_local_vertices_; ++i) {
    neighbors_of_vertices_offsets_[i] = j;
    j += ij;
  }

  j = 0;

  for (i = 0; i <= number_of_hyperedges_; ++i) {
    hyperedge_vertices_in_part_offsets_[i] = j;
    j += ij;
  }
}

void k_way_greedy_refiner::reset_data_structures() {
  to_non_local_vertices_.destroy();
  free_memory();
}

void k_way_greedy_refiner::set_partitioning_structures(int pNo, MPI_Comm comm) {
#ifdef DEBUG_REFINER
  assert(pNo >= 0 && pNo < number_of_partitions_);
#endif

  int i;
  int j;
  int ij;

  int neighVertOffset;
  int hEdgeVertOffset;
  int vertOffset;
  int endOffset;
  int numSpanned;
  int vertex;
  int vPart;
  int v;
  int nonLocIndex;

#ifdef DEBUG_REFINER
  assert(numParts > 1);
  assert(numNeighParts.getLength() > 0);
  assert(neighboursOfVOffsets.getLength() > 0);
  assert(neighboursOfV.getLength() > 0);
  assert(hEdgeVinPartOffsets.getLength() > 0);
#endif

  dynamic_array<int> locPartWts(number_of_parts_);

  // ###
  // initialise current partition vector
  // ###

  current_partition_vector_ = &partition_vector_[partition_vector_offsets_[pNo]];
  if (number_of_non_local_vertices_ == 0)
    current_non_local_partition_vector_ = nullptr;
  else
    current_non_local_partition_vector_ = &part_indices_[index_into_part_indices_[pNo]];
  current_partition_number_ = pNo;

  for (i = 0; i < number_of_parts_; ++i)
    locPartWts[i] = 0;

  // ###
  // initialise other data_ structures
  // first reset to zeros...
  // ###

  j = number_of_local_vertices_;

  for (i = 0; i < j; ++i) {
    locPartWts[current_partition_vector_[i]] += vertex_weights_[i];
    number_of_neighbor_parts_[i] = 0;
  }

  j = neighbors_of_vertices_offsets_[number_of_local_vertices_];

  for (i = 0; i < j; ++i)
    neighbors_of_vertices_[i] = 0;

  j = number_of_hyperedges_;

  for (i = 0; i < j; ++i)
    number_of_parts_spanned_[i] = 0;

  j = hyperedge_vertices_in_part_offsets_[number_of_hyperedges_];

  for (i = 0; i < j; ++i)
    hyperedge_vertices_in_part_[i] = 0;

  MPI_Allreduce(locPartWts.data(), part_weights_.data(), number_of_parts_,
                MPI_INT, MPI_SUM, comm);

#ifdef DEBUG_REFINER
  for (i = 0; i < numParts; ++i)
    assert(partWeights[i] <= maxPartWt);
#endif

  for (i = 0; i < number_of_hyperedges_; ++i) {
    endOffset = hyperedge_offsets_[i + 1];
    hEdgeVertOffset = hyperedge_vertices_in_part_offsets_[i];
    numSpanned = 0;

    for (j = hyperedge_offsets_[i]; j < endOffset; ++j) {
      ij = local_pin_list_[j];

      if (ij >= minimum_vertex_index_ && ij < maximum_vertex_index_) {
        vPart = current_partition_vector_[ij - minimum_vertex_index_];
      } else {
        nonLocIndex = to_non_local_vertices_.get(ij);
#ifdef DEBUG_REFINER
        assert(nonLocIndex >= 0 && nonLocIndex < numNonLocVerts);
#endif
        vPart = current_non_local_partition_vector_[nonLocIndex];
      }
#ifdef DEBUGU_REFINER
      assert(vPart >= 0 && vPart < numParts);
#endif
      ij = hEdgeVertOffset + vPart;

      if (hyperedge_vertices_in_part_[ij] == 0) {
        spanned_parts_[numSpanned++] = vPart;
      }

      ++hyperedge_vertices_in_part_[ij];
    }

#ifdef DEBUG_REFINER
    int hEdgeLen = 0;
    for (int ijk = 0; ijk < numParts; ++ijk)
      hEdgeLen += hEdgeVinPart[hEdgeVertOffset + ijk];
    assert(hEdgeLen == hEdgeOffset[i + 1] - hEdgeOffset[i]);
#endif

    number_of_parts_spanned_[i] = numSpanned;

    // ###
    // update the vertex neighbour structs
    // ###

    for (j = hyperedge_offsets_[i]; j < endOffset; ++j) {
      vertex = local_pin_list_[j];

      // ###
      // only consider if vertex a local vertex
      // ###

      if (vertex >= minimum_vertex_index_ && vertex < maximum_vertex_index_) {
        v = vertex - minimum_vertex_index_;
        vertOffset = neighbors_of_vertices_offsets_[v];

        for (ij = 0; ij < numSpanned; ++ij) {
          neighVertOffset = vertOffset + spanned_parts_[ij];

          if (neighbors_of_vertices_[neighVertOffset] == 0)
            ++number_of_neighbor_parts_[v];

          neighbors_of_vertices_[neighVertOffset] = 1;
        }
      }
    }
  }

#ifdef DEBUG_REFINER
  nonLocVertCheck();
  sanityHedgeCheck();
#endif
}

void k_way_greedy_refiner::refine(parallel::hypergraph &h,
                                           MPI_Comm comm) {
  initialize_data_structures(h, comm);

  int i;

  int gain;
  int lastGain = 0;
  int totalGain;
  int newCutsize;
  int numPasses;

  for (i = 0; i < number_of_partitions_; ++i) {
    set_partitioning_structures(i, comm);

    totalGain = 0;
    numPasses = 0;

    progress("\t[%i] ", i);

    do {
      newCutsize = greedy_k_way_refinement(h, i, comm);
      gain = partition_cuts_[i] - newCutsize;

      progress("%i ", gain);
      if (gain < 0) {
        undo_pass_moves();
      } else {
        partition_cuts_[i] = newCutsize;
        totalGain += gain;
        ++numPasses;
      }

      if (early_exit_ && numPasses > 1 && gain < lastGain)
        break;

      lastGain = gain;
    } while (gain > 0);
    progress("| %i %i %i\n", numPasses, totalGain, partition_cuts_[i]);
  }

  reset_data_structures();
}

int k_way_greedy_refiner::greedy_k_way_refinement(parallel::hypergraph &h,
                                                  int pNo, MPI_Comm comm) {
  int i;

  total_number_of_vertices_moved_ = 0;
  locked_.unset();

  for (i = 0; i < 2; ++i) {
    if (rank_ == ROOT_PROC) {
      // ###
      // init movement sets
      // ###
      movement_sets_->initialize_part_weights(part_weights_.data(),
                                              number_of_parts_);
    }

    greedy_pass(i, comm);
    manage_balance_constraint(comm);
    update_vertex_move_info(comm);
  }

  if (percentile_ == 100)
    return (compute_cutsize(comm));
  else {
    return (h.calculate_cut_size(number_of_parts_, pNo, comm));
  }
}

int k_way_greedy_refiner::greedy_pass(int lowToHigh, MPI_Comm comm) {
  int i;
  int j;
  int ij;
  int v;
  int sP;
  int gain;
  int prod;
  int vGain;
  int posGain;
  int hEdge;
  int bestMove;
  int randomNum;
  int vertexWt;
  int vertOffset;
  int vNeighOffset;
  int hEdgeOff;
  int neighOfVOffset;
  int numNonPos = 0;
  int limNonPosMoves = static_cast<int>(
      ceil(limit_ * static_cast<double>(number_of_local_vertices_)));

  double currImbalance = 0;
  double bestImbalance;
  double posImbalance;

  for (i = 0; i < number_of_parts_; ++i) {
    prod = number_of_parts_ * i;

    for (j = 0; j < number_of_parts_; ++j) {
      if (i != j) {
        ij = prod + j;

        number_of_vertices_moved_[ij] = 0;
        v = index_into_move_set_[ij];

        move_set_data_[v] = 0;
        move_set_data_[v + 1] = 0;
      }
    }
  }



  for (i = 0; i < number_of_local_vertices_; ++i) {
    vertices_[i] = i;
  }

  for (i = 0; i < number_of_parts_; ++i) {
    currImbalance += std::fabs(part_weights_[i] - average_part_weight_);
  }

  i = number_of_local_vertices_;
  gain = 0;

  do {
    randomNum = parkway::utility::random(0, i);
    v = vertices_[randomNum];

    if (!locked_(v)) {
      vertexWt = vertex_weights_[v];
      sP = current_partition_vector_[v];
      vGain = 0;
      bestImbalance = currImbalance;
      bestMove = -1;

      if (number_of_neighbor_parts_[v] > 1) {
        vNeighOffset = neighbors_of_vertices_offsets_[v];

        for (j = 0; j < number_of_parts_; ++j) {
          if (neighbors_of_vertices_[vNeighOffset + j] > 0 &&
              ((lowToHigh && j > sP) || (!lowToHigh && j < sP))) {
            if (part_weights_[j] + vertexWt <= maximum_part_weight_) {
              posGain = 0;
              vertOffset = vertex_to_hyperedges_offset_[v + 1];

              for (ij = vertex_to_hyperedges_offset_[v]; ij < vertOffset; ++ij) {
                hEdge = vertex_to_hyperedges_[ij];
                hEdgeOff = hyperedge_vertices_in_part_offsets_[hEdge];

                if (hyperedge_vertices_in_part_[hEdgeOff + sP] == 1) {
                  posGain += hyperedge_weights_[hEdge];
                }

                if (hyperedge_vertices_in_part_[hEdgeOff + j] == 0) {
                  posGain -= hyperedge_weights_[hEdge];
                }
              }

              posImbalance = currImbalance + std::fabs(part_weights_[sP] - (vertexWt + average_part_weight_));
              posImbalance += std::fabs((part_weights_[j] + vertexWt) - average_part_weight_);
              posImbalance -= std::fabs(part_weights_[sP] - average_part_weight_);
              posImbalance -= std::fabs(part_weights_[j] - average_part_weight_);

              if ((posGain > vGain) ||
                  (posGain == vGain && posImbalance < bestImbalance)) {
                vGain = posGain;
                bestMove = j;
                bestImbalance = posImbalance;
              }
            }
          }
        }


        if (bestMove != -1) {
          vertOffset = vertex_to_hyperedges_offset_[v + 1];
          neighOfVOffset = neighbors_of_vertices_offsets_[v];

          // ###
          // update the moved vertices' stats
          // ###

          if (neighbors_of_vertices_[neighOfVOffset + bestMove] == 0) {
            ++number_of_neighbor_parts_[v];
          }

          neighbors_of_vertices_[neighOfVOffset + bestMove] = 1;
          neighbors_of_vertices_[neighOfVOffset + sP] = 0;

          for (j = vertex_to_hyperedges_offset_[v]; j < vertOffset; ++j) {
            // ###
            // update the hyperedge stats: (vInPart etc.)
            // ###

            hEdge = vertex_to_hyperedges_[j];
            hEdgeOff = hyperedge_vertices_in_part_offsets_[hEdge];
            if (hyperedge_vertices_in_part_[hEdgeOff + bestMove] == 0) {
              ++number_of_parts_spanned_[hEdge];
            }

            --hyperedge_vertices_in_part_[hEdgeOff + sP];
            ++hyperedge_vertices_in_part_[hEdgeOff + bestMove];

            if (hyperedge_vertices_in_part_[hEdgeOff + sP] > 0) {
              neighbors_of_vertices_[neighOfVOffset + sP] = 1;
            } else {
              --number_of_parts_spanned_[hEdge];
            }
          }

          if (neighbors_of_vertices_[neighOfVOffset + sP] == 0) {
            --number_of_neighbor_parts_[v];
          }

          // ###
          // update the adj vertices stats:
          // (num neighbours in part etc.)
          // ###

          update_adjacent_vertex_status(v, sP, bestMove);

          // ###
          // update other structs
          // ###

          locked_.set(v);
          current_partition_vector_[v] = bestMove;
          part_weights_[sP] -= vertexWt;
          part_weights_[bestMove] += vertexWt;
          currImbalance = bestImbalance;

          // ###
          // update the gain...
          // ###

          gain += vGain;
          numNonPos = vGain <= 0 ? numNonPos + 1 : 0;

          // ###
          // update the movement set structures
          // ###

          ij = sP * number_of_parts_ + bestMove;

          move_sets_[ij]->at(number_of_vertices_moved_[ij]++) = v + minimum_vertex_index_;
          move_set_data_[index_into_move_set_[ij]] += vGain;
          move_set_data_[index_into_move_set_[ij] + 1] += vertexWt;

          if (limit_ < 1.0 && numNonPos > limNonPosMoves) {
            break;
          }
        }
      }
    }

    std::swap(vertices_[randomNum], vertices_[--i]);
  } while (i > 0);
  return gain;
}

int k_way_greedy_refiner::compute_cutsize(MPI_Comm comm) {
  int i;
  int ij;
  int locCut = 0;

  int totalCut;

  for (i = 0; i < number_of_allocated_hyperedges_; ++i) {
    ij = allocated_hyperedges_[i];
#ifdef DEBUG_REFINER
    assert(numPartsSpanned[ij] > 0 && numPartsSpanned[ij] <= numParts);
    assert(hEdgeWeight[ij] > 0);
#endif
    locCut += ((number_of_parts_spanned_[ij] - 1) * hyperedge_weights_[ij]);
  }

  MPI_Allreduce(&locCut, &totalCut, 1, MPI_INT, MPI_SUM, comm);

  return totalCut;
}

void k_way_greedy_refiner::manage_balance_constraint(MPI_Comm comm) {
  int i;
  int j;
  int ij;

  int prod;
  int arrayLen;
  int numToSend;
  int totToRecv;
  int indexIntoMoveSets;

  ds::dynamic_array<int> movesLengths;

  dynamic_array<dynamic_array<int> *> moves;

  numToSend = 0;
  for (i = 0; i < number_of_parts_; ++i) {
    prod = i * number_of_parts_;

    for (j = 0; j < number_of_parts_; ++j) {
      if (i != j) {
        ij = index_into_move_set_[prod + j];

        if (move_set_data_[ij + 1] > 0) {
          send_array_[numToSend++] = i;
          send_array_[numToSend++] = j;
          send_array_[numToSend++] = move_set_data_[ij];
          send_array_[numToSend++] = move_set_data_[ij + 1];
        } else {
        }
      }
    }
  }

  MPI_Gather(&numToSend, 1, MPI_INT, receive_lens_.data(), 1, MPI_INT,
             ROOT_PROC, comm);

  if (rank_ == ROOT_PROC) {
    ij = 0;

    for (i = 0; i < processors_; ++i) {
      receive_displs_[i] = ij;
      ij += receive_lens_[i];
    }

    receive_array_.resize(ij);
    totToRecv = ij;
  }

  MPI_Gatherv(send_array_.data(), numToSend, MPI_INT, receive_array_.data(),
              receive_lens_.data(), receive_displs_.data(), MPI_INT, ROOT_PROC,
              comm);

  if (rank_ == ROOT_PROC) {
    for (i = 0; i < processors_; ++i) {
      ij = receive_displs_[i];

      if (receive_lens_[i] > 0) {
        movement_sets_->complete_processor_sets(i, receive_lens_[i],
                                              &(receive_array_[ij]));
      }
    }

    movement_sets_->compute_restoring_array();
  }

  moves = movement_sets_->restoring_moves();
  movesLengths = movement_sets_->restoring_move_lens();

  MPI_Scatter(movesLengths.data(), 1, MPI_INT, &totToRecv, 1, MPI_INT, ROOT_PROC,
              comm);

  if (rank_ == ROOT_PROC) {
    ij = 0;

    for (i = 0; i < processors_; ++i) {
      arrayLen = movesLengths[i];
      send_displs_[i] = ij;

      for (j = 0; j < arrayLen; ++j)
        send_array_[ij++] = moves[i]->at(j);;
    }
  }

  receive_array_.resize(totToRecv);

  MPI_Scatterv(send_array_.data(), movesLengths.data(), send_displs_.data(),
               MPI_INT, receive_array_.data(), totToRecv, MPI_INT, ROOT_PROC,
               comm);

  // ###
  // now everyone has been told which sets of moves to
  // take back in order to maintain balance condition
  // take back local moves as directed by root processor

  // note that can decrease the size of the take-back data_
  // since we are not making use of the weight of the set
  // to be taken back
  // ###

  ij = 0;

  while (ij < totToRecv) {
    i = receive_array_[ij];
    j = receive_array_[ij + 1];

    indexIntoMoveSets = i * number_of_parts_ + j;
    undo_move(indexIntoMoveSets, i, j);

    ij += 2;
  }

  if (rank_ == ROOT_PROC) {
    for (ij = 0; ij < number_of_parts_; ++ij) {
      part_weights_[ij] = movement_sets_->part_weights_array()[ij];
    }
  }

  MPI_Bcast(part_weights_.data(), number_of_parts_, MPI_INT, ROOT_PROC, comm);
}

void k_way_greedy_refiner::undo_pass_moves() {
  int i;
  int v;
  int vPart;

  for (i = 0; i < total_number_of_vertices_moved_; i += 2) {
    v = moved_vertices_[i];
    vPart = moved_vertices_[i + 1];
    current_partition_vector_[v - minimum_vertex_index_] = vPart;
  }
}

void k_way_greedy_refiner::update_vertex_move_info(MPI_Comm comm) {
#ifdef DEBUG_REFINER
  sanityHedgeCheck();
#endif

  int i;
  int j;
  int ij;
  int ijk;

  int index;
  int prod;
  int vert;
  int totToSend;
  int totToRecv;
  int minIndex;
  int v;
  int vertexPart;
  int newVertexPart;
  int locVertIndex;
  int hEdge;
  int endOffset;
  int neighOfVOffset;
  int othVOffset;
  int othHedge;
  int hEdgeOff;
  int numVerticesSeen;
  int nonLocIdx;
  int endNonLocOffset;

  totToSend = 0;
  for (i = 0; i < number_of_parts_; ++i) {
    prod = i * number_of_parts_;

    for (j = 0; j < number_of_parts_; ++j) {
      if (i != j) {
        index = prod + j;
        endOffset = number_of_vertices_moved_[index];
        for (ij = 0; ij < endOffset; ++ij) {
          v = move_sets_[index]->at(ij);
          send_array_[totToSend++] = v;
          send_array_[totToSend++] = j;

          moved_vertices_[total_number_of_vertices_moved_++] = v;
          moved_vertices_[total_number_of_vertices_moved_++] = i;
        }
      }
    }
  }

  MPI_Allgather(&totToSend, 1, MPI_INT, receive_lens_.data(), 1, MPI_INT, comm);

  ij = 0;
  for (i = 0; i < processors_; ++i) {
    receive_displs_[i] = ij;
    ij += receive_lens_[i];
  }

  receive_array_.resize(ij);
  totToRecv = ij;

  MPI_Allgatherv(send_array_.data(), totToSend, MPI_INT,
                 receive_array_.data(), receive_lens_.data(),
                 receive_displs_.data(), MPI_INT, comm);

  // ###
  // now go through the moved vertices and update local structures
  // ###

  minIndex = receive_displs_[rank_];

#ifdef DEBUG_REFINER
  assert(And(totToRecv, 0x1) == 0);
#endif

  i = 0;

  while (i < minIndex) {
    v = receive_array_[i];

#ifdef DEBUG_REFINER
    assert(v < minVertexIndex || v >= maxVertexIndex);
#endif
    nonLocIdx = to_non_local_vertices_.get_careful(v);

    if (nonLocIdx >= 0) {
#ifdef DEBUG_REFINER
      assert(nonLocIdx < numNonLocVerts);
#endif
      vertexPart = current_non_local_partition_vector_[nonLocIdx];
      newVertexPart = receive_array_[i + 1];
      endNonLocOffset = non_local_vertices_to_hyperedges_offsets_[nonLocIdx + 1];

#ifdef DEBUG_REFINER
      assert(newVertexPart >= 0 && newVertexPart < numParts);
      assert(vertexPart >= 0 && vertexPart < numParts);
      assert(vertexPart != newVertexPart);
#endif

      // ###
      // first update the hyperedge vInPart structure
      // ###

      for (j = non_local_vertices_to_hyperedges_offsets_[nonLocIdx]; j < endNonLocOffset; ++j) {
        hEdge = non_local_vertices_to_hyperedges_[j];

#ifdef DEBUG_REFINER
        assert(hEdge >= 0 && hEdge < numHedges);
#endif
        hEdgeOff = hyperedge_vertices_in_part_offsets_[hEdge];

#ifdef DEBUG_REFINER
        int hEdgeLen = 0;
        for (int ijkl = 0; ijkl < numParts; ++ijkl)
          hEdgeLen += hEdgeVinPart[hEdgeOff + ijkl];
        assert(hEdgeLen == hEdgeOffset[hEdge + 1] - hEdgeOffset[hEdge]);
        assert(hEdgeVinPart[hEdgeOff + vertexPart] > 0);
#endif
        if (hyperedge_vertices_in_part_[hEdgeOff + newVertexPart] == 0) {
          ++number_of_parts_spanned_[hEdge];
        }

        --hyperedge_vertices_in_part_[hEdgeOff + vertexPart];
        ++hyperedge_vertices_in_part_[hEdgeOff + newVertexPart];

#ifdef DEBUG_REFINER
        assert(hEdgeVinPart[hEdgeOff + newVertexPart] > 0);
        assert(hEdgeVinPart[hEdgeOff + vertexPart] >= 0);
#endif

        if (hyperedge_vertices_in_part_[hEdgeOff + vertexPart] == 0)
          --number_of_parts_spanned_[hEdge];
      }

      // ###
      // update the data_ for local adjacent vertices
      // ###

      numVerticesSeen = 0;

      for (j = non_local_vertices_to_hyperedges_offsets_[nonLocIdx]; j < endNonLocOffset; ++j) {
        hEdge = non_local_vertices_to_hyperedges_[j];
#ifdef DEBUG_REFINER
        assert(hEdge >= 0 && hEdge < numHedges);
#endif
        hEdgeOff = hyperedge_offsets_[hEdge + 1];

        for (ij = hyperedge_offsets_[hEdge]; ij < hEdgeOff; ++ij) {
          vert = local_pin_list_[ij];

          // ###
          // if the adjacent vertex is a local vertex
          // ###

          if (vert >= minimum_vertex_index_ && vert < maximum_vertex_index_) {
            locVertIndex = vert - minimum_vertex_index_;

            if (vertices_seen_(locVertIndex) == 0) {
              neighOfVOffset = neighbors_of_vertices_offsets_[locVertIndex];

              if (neighbors_of_vertices_[neighOfVOffset + newVertexPart] == 0) {
                ++number_of_neighbor_parts_[locVertIndex];
              }

              neighbors_of_vertices_[neighOfVOffset + newVertexPart] = 1;

              if (current_partition_vector_[locVertIndex] != vertexPart) {
                neighbors_of_vertices_[neighOfVOffset + vertexPart] = 0;

                // sorting out this now...go thro each hedge
                // of neighbouring vertex and check if its
                // connected to part bestDP to determine
                // whether to set to 0 or 1!!!

                othVOffset = vertex_to_hyperedges_offset_[locVertIndex + 1];

                for (ijk = vertex_to_hyperedges_offset_[locVertIndex]; ijk < othVOffset;
                     ++ijk) {
                  othHedge = vertex_to_hyperedges_[ijk];

                  if (hyperedge_vertices_in_part_[
                          hyperedge_vertices_in_part_offsets_[othHedge] + vertexPart] >
                      0) {
                    neighbors_of_vertices_[neighOfVOffset + vertexPart] = 1;
                    break;
                  }
                }

                if (neighbors_of_vertices_[neighOfVOffset + vertexPart] == 0) {
                  --number_of_neighbor_parts_[locVertIndex];
                }
              }

              vertices_seen_.set(locVertIndex);
              seen_vertices_[numVerticesSeen++] = locVertIndex;
            }
          }

          // ###
          // if the adjacent vertex is not a local vertex
          // we don't update anything
          // ###
        }
      }

      current_non_local_partition_vector_[nonLocIdx] = newVertexPart;

      // ###
      // restore the 'seen' vertices structure
      // ###

      for (j = 0; j < numVerticesSeen; ++j) {
        vertices_seen_.unset(seen_vertices_[j]);
      }

      // ###
      // end updating structures after remote vertex move
      // ###
    }

    i += 2;
  }

  i += totToSend;

  while (i < totToRecv) {
    v = receive_array_[i];

#ifdef DEBUG_REFINER
    assert(v < minVertexIndex || v >= maxVertexIndex);
#endif

    nonLocIdx = to_non_local_vertices_.get_careful(v);

    if (nonLocIdx >= 0) {
#ifdef DEBUG_REFINER
      assert(nonLocIdx < numNonLocVerts);
#endif
      vertexPart = current_non_local_partition_vector_[nonLocIdx];
      newVertexPart = receive_array_[i + 1];
      endNonLocOffset = non_local_vertices_to_hyperedges_offsets_[nonLocIdx + 1];

#ifdef DEBUG_REFINER
      assert(newVertexPart >= 0 && newVertexPart < numParts);
      assert(vertexPart >= 0 && vertexPart < numParts);
      assert(vertexPart != newVertexPart);
#endif

      // ###
      // first update the hyperedge vInPart structure
      // ###

      for (j = non_local_vertices_to_hyperedges_offsets_[nonLocIdx]; j < endNonLocOffset; ++j) {
        hEdge = non_local_vertices_to_hyperedges_[j];

#ifdef DEBUG_REFINER
        assert(hEdge >= 0 && hEdge < numHedges);
#endif
        hEdgeOff = hyperedge_vertices_in_part_offsets_[hEdge];
#ifdef DEBUG_REFINER
        int hEdgeLen = 0;
        for (int ijkl = 0; ijkl < numParts; ++ijkl)
          hEdgeLen += hEdgeVinPart[hEdgeOff + ijkl];
        assert(hEdgeLen == hEdgeOffset[hEdge + 1] - hEdgeOffset[hEdge]);
        assert(hEdgeVinPart[hEdgeOff + vertexPart] > 0);
#endif
        if (hyperedge_vertices_in_part_[hEdgeOff + newVertexPart] == 0) {
          ++number_of_parts_spanned_[hEdge];
        }

        --hyperedge_vertices_in_part_[hEdgeOff + vertexPart];
        ++hyperedge_vertices_in_part_[hEdgeOff + newVertexPart];

#ifdef DEBUG_REFINER
        assert(hEdgeVinPart[hEdgeOff + newVertexPart] > 0);
        assert(hEdgeVinPart[hEdgeOff + vertexPart] >= 0);
#endif
        if (hyperedge_vertices_in_part_[hEdgeOff + vertexPart] == 0)
          --number_of_parts_spanned_[hEdge];
      }

      // ###
      // update the data_ for local adjacent vertices
      // ###

      numVerticesSeen = 0;

      for (j = non_local_vertices_to_hyperedges_offsets_[nonLocIdx]; j < endNonLocOffset; ++j) {
        hEdge = non_local_vertices_to_hyperedges_[j];
#ifdef DEBUG_REFINER
        assert(hEdge >= 0 && hEdge < numHedges);
#endif
        hEdgeOff = hyperedge_offsets_[hEdge + 1];

        for (ij = hyperedge_offsets_[hEdge]; ij < hEdgeOff; ++ij) {
          vert = local_pin_list_[ij];

          // ###
          // if the adjacent vertex is a local vertex
          // ###

          if (vert >= minimum_vertex_index_ && vert < maximum_vertex_index_) {
            locVertIndex = vert - minimum_vertex_index_;

            if (vertices_seen_(locVertIndex) == 0) {
              neighOfVOffset = neighbors_of_vertices_offsets_[locVertIndex];

              if (neighbors_of_vertices_[neighOfVOffset + newVertexPart] == 0) {
                ++number_of_neighbor_parts_[locVertIndex];
              }

              neighbors_of_vertices_[neighOfVOffset + newVertexPart] = 1;

              if (current_partition_vector_[locVertIndex] != vertexPart) {
                neighbors_of_vertices_[neighOfVOffset + vertexPart] = 0;

                // sorting out this now...go thro each hedge
                // of neighbouring vertex and check if its
                // connected to part bestDP to determine
                // whether to set to 0 or 1!!!

                othVOffset = vertex_to_hyperedges_offset_[locVertIndex + 1];

                for (ijk = vertex_to_hyperedges_offset_[locVertIndex]; ijk < othVOffset;
                     ++ijk) {
                  othHedge = vertex_to_hyperedges_[ijk];

                  if (hyperedge_vertices_in_part_[
                          hyperedge_vertices_in_part_offsets_[othHedge] + vertexPart] >
                      0) {
                    neighbors_of_vertices_[neighOfVOffset + vertexPart] = 1;
                    break;
                  }
                }

                if (neighbors_of_vertices_[neighOfVOffset + vertexPart] == 0) {
                  --number_of_neighbor_parts_[locVertIndex];
                }
              }

              vertices_seen_.set(locVertIndex);
              seen_vertices_[numVerticesSeen++] = locVertIndex;
            }
          }

          // ###
          // if the adjacent vertex is not a local vertex
          // we don't update anything
          // ###
        }
      }

      current_non_local_partition_vector_[nonLocIdx] = newVertexPart;

      // ###
      // restore the 'seen' vertices structure
      // ###

      for (j = 0; j < numVerticesSeen; ++j)
        vertices_seen_.unset(seen_vertices_[j]);

      // ###
      // end updating structures after remote vertex move
      // ###
    }

    i += 2;
  }

#ifdef DEBUG_REFINER
  sanityHedgeCheck();
#endif
}

void k_way_greedy_refiner::update_adjacent_vertex_status(int v, int sP,
                                                         int bestMove) {
  int i;
  int j;
  int ij;

  int numVerticesSeen;
  int vertOffset;
  int vert;
  int locVertIndex;
  int hEdge;
  int hEdgeOff;
  int neighOfVOffset;
  int othVOffset;
  int othHedge;

  vertOffset = vertex_to_hyperedges_offset_[v + 1];
  numVerticesSeen = 0;

  for (i = vertex_to_hyperedges_offset_[v]; i < vertOffset; ++i) {
    hEdge = vertex_to_hyperedges_[i];
    hEdgeOff = hyperedge_offsets_[hEdge + 1];

    for (j = hyperedge_offsets_[hEdge]; j < hEdgeOff; ++j) {
      vert = local_pin_list_[j];

      // ###
      // if the adjacent vertex is a local vertex
      // ###

      if (vert >= minimum_vertex_index_ && vert < maximum_vertex_index_) {
        locVertIndex = vert - minimum_vertex_index_;

        if (locVertIndex != v && vertices_seen_(locVertIndex) == 0) {
          neighOfVOffset = neighbors_of_vertices_offsets_[locVertIndex];

          if (neighbors_of_vertices_[neighOfVOffset + bestMove] == 0) {
            ++number_of_neighbor_parts_[locVertIndex];
          }

          neighbors_of_vertices_[neighOfVOffset + bestMove] = 1;

          if (current_partition_vector_[locVertIndex] != sP) {
            neighbors_of_vertices_[neighOfVOffset + sP] = 0;

            // sorting out this now...go thro each hedge
            // of neighbouring vertex and check if its
            // connected to part bestDP to determine
            // whether to set to 0 or 1!!!

            othVOffset = vertex_to_hyperedges_offset_[locVertIndex + 1];

            for (ij = vertex_to_hyperedges_offset_[locVertIndex]; ij < othVOffset; ++ij) {
              othHedge = vertex_to_hyperedges_[ij];

              if (hyperedge_vertices_in_part_[
                      hyperedge_vertices_in_part_offsets_[othHedge] + sP] > 0) {
                neighbors_of_vertices_[neighOfVOffset + sP] = 1;
                break;
              }
            }

            if (neighbors_of_vertices_[neighOfVOffset + sP] == 0) {
              --number_of_neighbor_parts_[locVertIndex];
            }
          }

          vertices_seen_.set(locVertIndex);
          seen_vertices_[numVerticesSeen++] = locVertIndex;
        }
      }

      // ###
      // if the adjacent vertex is not a local vertex
      // we don't update anything
      // ###
    }
  }

  // ###
  // restore the 'seen' vertices structure
  // ###

  for (i = 0; i < numVerticesSeen; ++i)
    vertices_seen_.unset(seen_vertices_[i]);
}

void k_way_greedy_refiner::undo_move(int indexIntoMoveSets, int from,
                                              int to) {
#ifdef DEBUG_REFINER
  sanityHedgeCheck();
#endif

  int v;
  int j;

  int i;
  int hEdgeOff;
  int neighOfVOffset;
  int vertOffset;
  int hEdge;

  int numVmoved = number_of_vertices_moved_[indexIntoMoveSets];
  int *movedArray = move_sets_[indexIntoMoveSets]->data();

  for (i = 0; i < numVmoved; ++i) {
    v = movedArray[i] - minimum_vertex_index_;

#ifdef DEBUG_REFINER
    assert(v >= 0 && v < numLocalVertices);
#endif

    // ###
    // make the move and update the vertex's structs
    // ###

    vertOffset = vertex_to_hyperedges_offset_[v + 1];
    neighOfVOffset = neighbors_of_vertices_offsets_[v];

    // ###
    // update the moved vertices' stats
    // ###

    if (neighbors_of_vertices_[neighOfVOffset + from] == 0) {
      ++number_of_neighbor_parts_[v];
    }

    neighbors_of_vertices_[neighOfVOffset + from] = 1;
    neighbors_of_vertices_[neighOfVOffset + to] = 0;

    for (j = vertex_to_hyperedges_offset_[v]; j < vertOffset; ++j) {
      // ###
      // update the hyperedge stats: (vInPart etc.)
      // ###

      hEdge = vertex_to_hyperedges_[j];
#ifdef DEBUG_REFINER
      assert(hEdge >= 0 && hEdge < numHedges);
#endif
      hEdgeOff = hyperedge_vertices_in_part_offsets_[hEdge];

#ifdef DEBUG_REFINER
      int hEdgeLen = 0;
      for (int ijk = 0; ijk < numParts; ++ijk)
        hEdgeLen += hEdgeVinPart[hEdgeOff + ijk];
      assert(hEdgeLen == hEdgeOffset[hEdge + 1] - hEdgeOffset[hEdge]);
#endif
      if (hyperedge_vertices_in_part_[hEdgeOff + from] == 0) {
        ++number_of_parts_spanned_[hEdge];
      }

      --hyperedge_vertices_in_part_[hEdgeOff + to];
      ++hyperedge_vertices_in_part_[hEdgeOff + from];

#ifdef DEBUG_REFINER
      assert(hEdgeVinPart[hEdgeOff + to] >= 0);
#endif

      if (hyperedge_vertices_in_part_[hEdgeOff + to] > 0) {
        neighbors_of_vertices_[neighOfVOffset + to] = 1;
      } else {
        --number_of_parts_spanned_[hEdge];
      }
    }

    if (neighbors_of_vertices_[neighOfVOffset + to] == 0) {
      --number_of_neighbor_parts_[v];
    }

    // ###
    // update the adj vertices stats:
    // (num neighbours in part etc.)
    // ###

    update_adjacent_vertex_status(v, to, from);

    // ###
    // update other structs
    // ###

    locked_.unset(v);
    current_partition_vector_[v] = from;
  }

#ifdef DEBUG_REFINER
  sanityHedgeCheck();
#endif

  number_of_vertices_moved_[indexIntoMoveSets] = 0;
}

void k_way_greedy_refiner::non_local_vertices_check() const {
  int i;
  int j;
  int ij;

  int vertPart;
  int endOffset;
  int h;

  for (i = 0; i < number_of_non_local_vertices_; ++i) {
    vertPart = current_non_local_partition_vector_[i];
    endOffset = non_local_vertices_to_hyperedges_offsets_[i + 1];

    for (j = non_local_vertices_to_hyperedges_offsets_[i]; j < endOffset; ++j) {
      assert(j < non_local_vertices_to_hyperedges_.capacity());
      h = non_local_vertices_to_hyperedges_[j];
      assert(h < hyperedge_vertices_in_part_offsets_.capacity());
      ij = hyperedge_vertices_in_part_offsets_[h];
      assert(hyperedge_vertices_in_part_[ij + vertPart] > 0);
    }
  }
}

void k_way_greedy_refiner::sanity_hyperedge_check() const {
  int i;
  int j;
  int ij;

  int hEdgeLen;
  int inParts;

  for (i = 0; i < number_of_hyperedges_; ++i) {
    hEdgeLen = hyperedge_offsets_[i + 1] - hyperedge_offsets_[i];
    inParts = 0;
    ij = hyperedge_vertices_in_part_offsets_[i + 1];

    for (j = hyperedge_vertices_in_part_offsets_[i]; j < ij; ++j) {
      inParts += hyperedge_vertices_in_part_[j];
      assert(hyperedge_vertices_in_part_[j] >= 0);
    }

    assert(inParts == hEdgeLen);
  }
}

}  // namespace parallel
}  // namespace parkway
