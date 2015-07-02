// ### Para2DModelCoarsener.cpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 19/04/2005: Last Modified
//
// ###
#include "coarseners/parallel/model_coarsener_2d.hpp"
#include "data_structures/map_to_pos_int.hpp"
#include "data_structures/bit_field.hpp"
#include "data_structures/internal/table_utils.hpp"
#include "utility/logging.hpp"

namespace parkway {
namespace parallel {

model_coarsener_2d::model_coarsener_2d(int rank, int nProcs, int nParts,
                                       int vertVisOrder, int matchReqOrder,
                                       int divByWt, int divByLen)
    : coarsener(rank, nProcs, nParts) {
  vertex_visit_order_ = vertVisOrder;
  match_request_visit_order_ = matchReqOrder;
  divide_by_cluster_weight_ = divByWt;
  divide_by_hyperedge_length_ = divByLen;
  limit_on_index_during_coarsening_ = 0;

  table_ = nullptr;
}

model_coarsener_2d::~model_coarsener_2d() {
}


void model_coarsener_2d::display_options() const {
  info("|--- PARA_C:\n"
       "|- 2DModel: r = %.2f min = %i vvo = ", reduction_ratio_,
       minimum_number_of_nodes_);
  print_visit_order(vertex_visit_order_);
  info(" mvo = ");
  print_visit_order(match_request_visit_order_);
  info(" divWt = %i divLen = %i threshold = %e\n|\n",
       divide_by_cluster_weight_, divide_by_hyperedge_length_, 1e6);
}

void model_coarsener_2d::build_auxiliary_structures(int numTotPins,
                                                    double aveVertDeg,
                                                    double aveHedgeSize) {
  // ###
  // build the ds::match_request_table
  // ###
  int i = static_cast<int>(ceil(static_cast<double>(numTotPins) / aveVertDeg));

  table_ = new ds::match_request_table(
      ds::internal::table_utils::table_size(i / processors_));
}

void model_coarsener_2d::release_memory() {
  hyperedge_weights_.reserve(0);
  hyperedge_offsets_.reserve(0);
  local_pin_list_.reserve(0);

  vertex_to_hyperedges_offset_.reserve(0);
  vertex_to_hyperedges_.reserve(0);
  allocated_hyperedges_.reserve(0);
  cluster_weights_.reserve(0);

    free_memory();
}

parallel::hypergraph *model_coarsener_2d::coarsen(
    parallel::hypergraph &h, MPI_Comm comm) {
  load(h, comm);

  if (number_of_vertices_ < minimum_number_of_nodes_ || h.dont_coarsen())
    return nullptr;

  if (number_of_vertices_ < 3000000) {
    return (parallel_first_choice_coarsen(h, comm));
  } else {
    return (parallel_hyperedge_coarsen(h, comm));
  }
}

parallel::hypergraph *model_coarsener_2d::parallel_first_choice_coarsen(
    parallel::hypergraph &h, MPI_Comm comm) {
#ifdef MEM_CHECK
  MPI_Barrier(comm);
  write_log(rank_, "[begin PFCC]: usage: %f", MemoryTracker::usage());
  Funct::printMemUse(rank_, "[begin PFCC]");
#endif

  int i;
  int j;

  int index = 0;
  int numNotMatched = number_of_local_vertices_;
  int maxLocWt = 0;
  int maxWt;
  int aveVertexWt = static_cast<int>(
      ceil(static_cast<double>(total_hypergraph_weight_) / number_of_vertices_));
  int bestMatch;
  int bestMatchWt = -1;
  int numVisited;
  int vPerProc = number_of_vertices_ / processors_;
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

  if (number_of_local_pins_ < number_of_vertices_ / 2)
    matchInfoLoc.create(number_of_local_pins_, 1);
  else
    matchInfoLoc.create(number_of_vertices_, 0);

  dynamic_array<int> neighVerts;
  dynamic_array<int> neighCluWts;
  dynamic_array<double> connectVals;
  dynamic_array<int> vertices(number_of_local_vertices_);

  permute_vertices_arrays(vertices, number_of_local_vertices_);

  if (parkway::utility::status::handler::progress_enabled()) {
    for (i = 0; i < number_of_local_vertices_; ++i) {
      if (vertex_weights_[i] > maxLocWt)
        maxLocWt = vertex_weights_[i];
    }

    MPI_Reduce(&maxLocWt, &maxWt, 1, MPI_INT, MPI_MAX, 0, comm);
    progress(" %i %i %i ", maximum_vertex_weight_, maxWt, aveVertexWt);
  }

  metric = static_cast<double>(number_of_local_vertices_) / reduction_ratio_;
  limit_on_index_during_coarsening_ =
      number_of_local_vertices_ - static_cast<int>(floor(metric - 1.0));
  cluster_index_ = 0;
  stop_coarsening_ = 0;

  // write_log(rank_, "numLocalVertices = %d, numLocalPins = %d, numHedges =
  // %d", numLocalVertices, numLocPins, numHedges);

  for (; index < number_of_local_vertices_; ++index) {
    if (index % 50000 == 0 && rank_ == 0)
      write_log(rank_, "considering local vertex index %d", index);

    if (match_vector_[vertices[index]] == -1) {
      vertex = vertices[index];
#ifdef DEBUG_COARSENER
      assert(vertex >= 0 && vertex < numLocalVertices);
#endif
      globalVertexIndex = vertex + minimum_vertex_index_;
      endOffset1 = vertex_to_hyperedges_offset_[vertex + 1];
      bestMatch = -1;
      maxMatchMetric = 0.0;
      numVisited = 0;

      for (i = vertex_to_hyperedges_offset_[vertex]; i < endOffset1; ++i) {
        hEdge = vertex_to_hyperedges_[i];
        endOffset2 = hyperedge_offsets_[hEdge + 1];
        hEdgeLen = endOffset2 - hyperedge_offsets_[hEdge];

        for (j = hyperedge_offsets_[hEdge]; j < endOffset2; ++j) {

          candidatV = local_pin_list_[j];
#ifdef DEBUG_COARSENER
          assert(candidatV >= 0 && candidatV < totalVertices);
#endif
          if (candidatV != globalVertexIndex) {
            /* now try not to check the weight before checking the vertex */

            neighbourLoc = matchInfoLoc.get_careful(candidatV);

            if (neighbourLoc >= 0) {
              if (divide_by_hyperedge_length_)
                connectVals[neighbourLoc] +=
                    (static_cast<double>(hyperedge_weights_[hEdge]) / (hEdgeLen - 1));
              else
                connectVals[neighbourLoc] +=
                    (static_cast<double>(hyperedge_weights_[hEdge]));
            } else {
              // ###
              // here compute the cluster weight
              // ###

              if (candidatV >= minimum_vertex_index_ && candidatV <
                                                        maximum_vertex_index_) {
                // ###
                // candidatV is a local vertex
                // ###

                if (match_vector_[candidatV - minimum_vertex_index_] == -1)
                  cluWeight =
                      vertex_weights_[vertex] + vertex_weights_[candidatV -
                                                minimum_vertex_index_];
                else if (match_vector_[candidatV - minimum_vertex_index_] >=
                         NON_LOCAL_MATCH) {
                  nonLocV =
                      match_vector_[candidatV - minimum_vertex_index_] - NON_LOCAL_MATCH;
                  cluWeight = vertex_weights_[vertex] +
                              table_->cluster_weight(nonLocV) + aveVertexWt;
                } else
                  cluWeight =
                      vertex_weights_[vertex] +
                      cluster_weights_[match_vector_[candidatV -
                                                 minimum_vertex_index_]];
              } else {
                // ###
                // candidatV is not a local vertex or it is a local
                // vertex matched with a non-local one
                // ###

                candVwt = table_->cluster_weight(candidatV);

                if (candVwt != -1)
                  cluWeight = vertex_weights_[vertex] + candVwt + aveVertexWt;
                else
                  cluWeight = vertex_weights_[vertex] + aveVertexWt;
              }

#ifdef DEBUG_COARSENER
              assert(numVisited >= 0);
#endif
              if (matchInfoLoc.insert(candidatV, numVisited)) {
                write_log(rank_, "numEntries %d",
                          matchInfoLoc.size());
                write_log(rank_, "using hash %d",
                          matchInfoLoc.use_hash());
                write_log(rank_, "numSlots %d", matchInfoLoc.capacity());
                write_log(rank_, "candidatV %d", candidatV);
                write_log(rank_, "neighbourLoc %d", neighbourLoc);
                assert(0);
              }

              neighVerts[numVisited] = candidatV;
              neighCluWts[numVisited] = cluWeight;

              if (divide_by_hyperedge_length_)
                connectVals[numVisited++] = static_cast<double>(hyperedge_weights_[hEdge]) / (hEdgeLen - 1);
              else
                connectVals[numVisited++] = static_cast<double>(hyperedge_weights_[hEdge]);
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

        if (pairWt <= maximum_vertex_weight_) {
          metric = connectVals[i];

          if (divide_by_cluster_weight_)
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

        match_vector_[vertex] = cluster_index_;
        cluster_weights_[cluster_index_++] = vertex_weights_[vertex];
        --numNotMatched;
      } else {
#ifdef DEBUG_COARSENER
        assert(bestMatch >= 0);
#endif

        if (bestMatch >= minimum_vertex_index_ && bestMatch <
                                                  maximum_vertex_index_) {
          // ###
          // best match is a local vertex
          // ###

          if (match_vector_[bestMatch - minimum_vertex_index_] == -1) {
            match_vector_[bestMatch - minimum_vertex_index_] = cluster_index_;
            match_vector_[vertex] = cluster_index_;
            cluster_weights_[cluster_index_++] = bestMatchWt;
            numNotMatched -= 2;
          } else {
            if (match_vector_[bestMatch - minimum_vertex_index_] >= NON_LOCAL_MATCH) {
              nonLocV =
                  match_vector_[bestMatch - minimum_vertex_index_] - NON_LOCAL_MATCH;
              table_->add_local(nonLocV, vertex + minimum_vertex_index_,
                               vertex_weights_[vertex],
                               std::min(nonLocV / vPerProc, processors_ - 1));
#ifdef DEBUG_COARSENER
              assert(std::min(nonLocV / vPerProc, processors_ - 1) != rank_);
#endif
              match_vector_[vertex] = NON_LOCAL_MATCH + nonLocV;
              --numNotMatched;
            } else {
              match_vector_[vertex] = match_vector_[bestMatch -
                                                minimum_vertex_index_];
              cluster_weights_[match_vector_[vertex]] +=
                  vertex_weights_[vertex]; // bestMatchWt;
              --numNotMatched;
            }
          }
        } else {
          // ###
          // best match is not a local vertex
          // ###

          table_->add_local(bestMatch, vertex + minimum_vertex_index_, vertex_weights_[vertex],
                           std::min(bestMatch / vPerProc, processors_ - 1));
#ifdef DEBUG_COARSENER
          assert(std::min(bestMatch / vPerProc, processors_ - 1) != rank_);
#endif
          match_vector_[vertex] = NON_LOCAL_MATCH + bestMatch;
          --numNotMatched;
        }
      }

      // ###
      // check the amount that the hypergraph has shrunk by
      // ###

      if (index > limit_on_index_during_coarsening_) {
        reducedBy = static_cast<double>(number_of_local_vertices_) /
                    (numNotMatched + cluster_index_ + table_->size());
        break;
      }
    }
  }

  // if(rank_ == 0)
  //  std::cout << "done local match" << std::endl;
  // MPI_Barrier(comm);

  matchInfoLoc.destroy();

  // ###
  // now carry over all the unmatched vertices as singletons
  // ###

  for (; index < number_of_local_vertices_; ++index) {
    vertex = vertices[index];

    if (match_vector_[vertex] == -1) {
      match_vector_[vertex] = cluster_index_;
      cluster_weights_[cluster_index_++] = vertex_weights_[vertex];
    }
  }

  // ###
  // now need to prepare and send out the matching requests
  // ###

  for (i = 0; i < 2; ++i) {
    set_request_arrays(i);
      send_from_data_out(comm); // actually sending requests
    set_reply_arrays(i, maximum_vertex_weight_);
      send_from_data_out(comm); // actually sending replies
    process_request_replies();
  }

  set_cluster_indices(comm);

  if (static_cast<double>(number_of_vertices_) / total_number_of_clusters_ <
      MIN_ALLOWED_REDUCTION_RATIO) {
    stop_coarsening_ = 1;
  }

  table_->clear();

  // ###
  // now construct the coarse hypergraph using the matching vector
  // ###
  return (contract_hyperedges(h, comm));
}

parallel::hypergraph *model_coarsener_2d::parallel_hyperedge_coarsen(
    parallel::hypergraph &h, MPI_Comm comm) {
#ifdef MEM_CHECK
  MPI_Barrier(comm);
  write_log(rank_, "[begin PFCC]: usage: %f", MemoryTracker::usage());
  Funct::printMemUse(rank_, "[begin PFCC]");
#endif

  int i;

  int index;
  int numNotMatched = number_of_local_vertices_;
  int maxLocWt = 0;
  int maxWt;
  int aveVertexWt = static_cast<int>(
      ceil(static_cast<double>(total_hypergraph_weight_) / number_of_vertices_));

  int vPerProc = number_of_vertices_ / processors_;

  int hEdge;
  int vertex;

  double reducedBy;

  if (parkway::utility::status::handler::progress_enabled()) {
    for (i = 0; i < number_of_local_vertices_; ++i) {
      if (vertex_weights_[i] > maxLocWt)
        maxLocWt = vertex_weights_[i];
    }

    MPI_Reduce(&maxLocWt, &maxWt, 1, MPI_INT, MPI_MAX, 0, comm);
    progress(" [PHEDGE] %i %i %i ", maximum_vertex_weight_, maxWt,
             aveVertexWt);
  }

  ds::bit_field matchedVertices(number_of_vertices_);
  matchedVertices.unset();

  dynamic_array<int> hEdges(number_of_hyperedges_);
  dynamic_array<int> hEdgeLens(number_of_hyperedges_);
  dynamic_array<int> unmatchedLocals;
  dynamic_array<int> unmatchedNonLocals;

  for (i = 0; i < number_of_hyperedges_; ++i) {
    hEdges[i] = i;
    hEdgeLens[i] = hyperedge_offsets_[i + 1] - hyperedge_offsets_[i];
  }

  /* sort hedges in increasing order of capacity_ */

  Funct::qsortByAnotherArray(0, number_of_hyperedges_ - 1, hEdges.data(),
                             hEdgeLens.data(), INC);

  /*
  */

  /* now sort the hyperedges of same capacity_ in
     decreasing order of weight */

  int leftBound = 0;
  int rightBound;
  int length;

  for (; leftBound < number_of_hyperedges_;) {
    length = hyperedge_weights_[leftBound];
    rightBound = leftBound + 1;

    for (; rightBound < number_of_hyperedges_ && hEdgeLens[rightBound] == length;
         ++rightBound)
      ;
    Funct::qsortByAnotherArray(leftBound, rightBound - 1, hEdges.data(),
                               hyperedge_weights_.data(), DEC);
    leftBound = rightBound;
  }

  cluster_index_ = 0;
  stop_coarsening_ = 0;

  for (index = 0; index < number_of_hyperedges_; ++index) {
    hEdge = hEdges[index];
    int numUnmatchedLocals = 0;
    int numUnmatchedNonLocals = 0;

    for (i = hyperedge_offsets_[hEdge]; i < hyperedge_offsets_[hEdge + 1]; ++i) {
      vertex = local_pin_list_[i];
#ifdef DEBUG_COARSENER
      assert(vertex >= 0 && vertex < totalVertices);
#endif
      if (vertex >= minimum_vertex_index_ && vertex < maximum_vertex_index_) {
        if (match_vector_[vertex - minimum_vertex_index_] == -1) {
          unmatchedLocals[numUnmatchedLocals++] = vertex - minimum_vertex_index_;
        }
      } else {
        if (matchedVertices(vertex) == 0) {
          unmatchedNonLocals[numUnmatchedNonLocals++] = vertex;
        }
      }
    }
    if (numUnmatchedLocals > 1 ||
        ((numUnmatchedLocals + numUnmatchedNonLocals) == hEdgeLens[hEdge])) {

      int totMatchWt = 0;

      for (i = 0; i < numUnmatchedLocals; ++i)
        totMatchWt += vertex_weights_[unmatchedLocals[i]];

      if (numUnmatchedLocals > 1) {
        if (totMatchWt < maximum_vertex_weight_) {
          for (i = 0; i < numUnmatchedLocals; ++i) {
            vertex = unmatchedLocals[i];
            match_vector_[vertex] = cluster_index_;
            cluster_weights_[cluster_index_] = totMatchWt;
          }

          ++cluster_index_;
          numNotMatched -= numUnmatchedLocals;
          reducedBy = static_cast<double>(number_of_local_vertices_) /
                      (numNotMatched + cluster_index_ + table_->size());

          if (reducedBy > reduction_ratio_)
            break;
        }
      } else {
        if (totMatchWt <
            (maximum_vertex_weight_ - (aveVertexWt * numUnmatchedNonLocals))) {
          int nonLocalVert = unmatchedNonLocals[0];

          for (i = 0; i < numUnmatchedLocals; ++i) {
            vertex = unmatchedLocals[i];
            table_->add_local(nonLocalVert, vertex + minimum_vertex_index_,
                             vertex_weights_[vertex],
                             std::min(nonLocalVert / vPerProc, processors_ - 1));
#ifdef DEBUG_COARSENER
            assert(std::min(nonLocalVert / vPerProc, processors_ - 1) != rank_);
#endif
            match_vector_[vertex] = NON_LOCAL_MATCH + nonLocalVert;
          }

          matchedVertices.set(nonLocalVert);

          numNotMatched -= numUnmatchedLocals;
          reducedBy = static_cast<double>(number_of_local_vertices_) /
                      (numNotMatched + cluster_index_ + table_->size());

          if (reducedBy > reduction_ratio_)
            break;
        }
      }
    }
  }

  for (index = 0; index < number_of_local_vertices_; ++index) {
    if (match_vector_[index] == -1) {
      match_vector_[index] = cluster_index_;
      cluster_weights_[cluster_index_++] = vertex_weights_[index];
    }
  }

  // ###
  // now need to prepare and send out the matching requests
  // ###

  for (i = 0; i < 2; ++i) {
    set_request_arrays(i);
      send_from_data_out(comm); // actually sending requests
    set_reply_arrays(i, maximum_vertex_weight_);
      send_from_data_out(comm); // actually sending replies
    process_request_replies();
  }

  set_cluster_indices(comm);

  if (static_cast<double>(number_of_vertices_) / total_number_of_clusters_ <
      MIN_ALLOWED_REDUCTION_RATIO) {
    stop_coarsening_ = 1;
  }

  table_->clear();

  // ###
  // now construct the coarse hypergraph using the matching vector
  // ###

  return (contract_hyperedges(h, comm));
}

void model_coarsener_2d::set_request_arrays(int highToLow) {
  int numRequests = table_->size();
  int nonLocVertex;
  int cluWt;

  int i;
  int procRank;

  ds::match_request_table::entry *entry;
  ds::dynamic_array<ds::match_request_table::entry *> entryArray = table_->get_entries();

  for (i = 0; i < processors_; ++i)
    send_lens_[i] = 0;

  for (i = 0; i < numRequests; ++i) {
    entry = entryArray[i];

#ifdef DEBUG_COARSENER
    assert(entry);
#endif

    nonLocVertex = entry->non_local_vertex();
    cluWt = entry->cluster_weight();
    procRank = entry->non_local_process();

#ifdef DEBUG_COARSENER
    assert(procRank < processors_);
#endif

    if ((cluWt >= 0) && ((highToLow && procRank < rank_) ||
                         (!highToLow && procRank > rank_))) {
      data_out_sets_[procRank][send_lens_[procRank]++] = nonLocVertex;
      data_out_sets_[procRank][send_lens_[procRank]++] = cluWt;
    }
  }
}

void model_coarsener_2d::set_reply_arrays(int highToLow, int maxVWt) {
  int j;
  int i;
  int l;

  dynamic_array<int> visitOrder;

  for (i = 0; i < processors_; ++i)
    send_lens_[i] = 0;

  int startOffset = 0;
  int vLocReq;
  int reqCluWt;
  int matchIndex;
  int visitOrderLen;

  for (i = 0; i < processors_; ++i) {
#ifdef DEBUG_COARSENER
    assert(And(receive_lens_[i], 0x1) == 0);
#endif

    if (match_request_visit_order_ == RANDOM_ORDER) {
      visitOrderLen = Shiftr(receive_lens_[i], 1);
      visitOrder.reserve(visitOrderLen);

      for (j = 0; j < visitOrderLen; ++j)
        visitOrder[j] = j;

      // ###
      // processing match requests in random order
      // ###

      Funct::randomPermutation(visitOrder.data(), visitOrderLen);

      for (l = 0; l < visitOrderLen; ++l) {

        j = Shiftl(visitOrder[l], 1);

        vLocReq = receive_array_[startOffset + j];
        reqCluWt = receive_array_[startOffset + j + 1];

        if (accept(vLocReq, reqCluWt, highToLow, maxVWt)) {
          // ###
          // cross-processor match accepted. inform vertices
          // of their cluster index and cluster weight
          // ###

          matchIndex = match_vector_[vLocReq - minimum_vertex_index_];
          data_out_sets_[i][send_lens_[i]++] = vLocReq;
          data_out_sets_[i][send_lens_[i]++] = matchIndex;
          data_out_sets_[i][send_lens_[i]++] = cluster_weights_[matchIndex];
        } else {
          // ###
          // cross-processor match rejected, inform vertices
          // that match rejected
          // ###

          data_out_sets_[i][send_lens_[i]++] = vLocReq;
          data_out_sets_[i][send_lens_[i]++] = NO_MATCH;
        }
      }
    } else {
      // ###
      // processing match requests as they arrive
      // ###

      visitOrderLen = receive_lens_[i];

      for (l = 0; l < visitOrderLen;) {
        vLocReq = receive_array_[startOffset + (l++)];
        reqCluWt = receive_array_[startOffset + (l++)];

        if (accept(vLocReq, reqCluWt, highToLow, maxVWt)) {
          // ###
          // cross-processor match accepted. inform vertices
          // of their cluster index and cluster weight
          // ###

          matchIndex = match_vector_[vLocReq - minimum_vertex_index_];
          data_out_sets_[i][send_lens_[i]++] = vLocReq;
          data_out_sets_[i][send_lens_[i]++] = matchIndex;
          data_out_sets_[i][send_lens_[i]++] = cluster_weights_[matchIndex];
        } else {
          // ###
          // cross-processor match rejected, inform vertices
          // that match rejected
          // ###

          data_out_sets_[i][send_lens_[i]++] = vLocReq;
          data_out_sets_[i][send_lens_[i]++] = NO_MATCH;
        }
      }
    }
    startOffset += receive_lens_[i];
  }
}

void model_coarsener_2d::process_request_replies() {
  int i;
  int j;
  int index;

  int startOffset = 0;
  int vNonLocReq;
  int cluWt;
  int matchIndex;
  int numLocals;

  ds::match_request_table::entry *entry;

  for (i = 0; i < processors_; ++i) {
    j = 0;
    while (j < receive_lens_[i]) {
      vNonLocReq = receive_array_[startOffset + (j++)];
      matchIndex = receive_array_[startOffset + (j++)];

      if (matchIndex != NO_MATCH) {
        // ###
        // match successful - set the clusterIndex
        // ###

        cluWt = receive_array_[startOffset + (j++)];
        entry = table_->get_entry(vNonLocReq);

#ifdef DEBUG_COARSENER
        assert(entry);
#endif

        entry->set_cluster_index(matchIndex);
        entry->set_cluster_weight(cluWt);
      } else {
        // ###
        // match not successful - match requesting
        // vertices into a cluster
        // ###

        entry = table_->get_entry(vNonLocReq);
        auto locals = entry->local_vertices_array();
        numLocals = entry->number_local();
        entry->set_cluster_index(MATCHED_LOCALLY);

        for (index = 0; index < numLocals; ++index)
          match_vector_[locals[index] - minimum_vertex_index_] = cluster_index_;

        cluster_weights_[cluster_index_++] = entry->cluster_weight();
      }
    }
    startOffset += receive_lens_[i];
  }
}

void model_coarsener_2d::set_cluster_indices(MPI_Comm comm) {
  dynamic_array<int> numClusters(processors_);
  dynamic_array<int> startIndex(processors_);

  MPI_Allgather(&cluster_index_, 1, MPI_INT, numClusters.data(), 1, MPI_INT,
                comm);

  int index = 0;
  int i;

  ds::match_request_table::entry *entry;
  ds::dynamic_array<ds::match_request_table::entry *> entryArray = table_->get_entries();

  int numLocals;
  int cluIndex;
  int numEntries = table_->size();

  i = 0;
  for (; index < processors_; ++index) {
    startIndex[index] = i;
    i += numClusters[index];
  }
  total_number_of_clusters_ = i;

  minimum_cluster_index_ = startIndex[rank_];

  for (index = 0; index < number_of_local_vertices_; ++index) {
#ifdef DEBUG_COARSENER
    assert(matchVector[index] != -1);
#endif
    if (match_vector_[index] < NON_LOCAL_MATCH)
      match_vector_[index] += minimum_cluster_index_;
  }

  // ###
  // now set the non-local matches' cluster indices
  // ###

  entryArray = table_->get_entries();

  for (index = 0; index < numEntries; ++index) {
    entry = entryArray[index];

#ifdef DEBUG_COARSENER
    assert(entry);
#endif

    cluIndex = entry->cluster_index();

#ifdef DEBUG_COARSENER
    assert(cluIndex >= 0);
#endif

    if (cluIndex != MATCHED_LOCALLY) {
      numLocals = entry->number_local();
      auto locals = entry->local_vertices_array();

#ifdef DEBUG_COARSENER
      assert(locals);
#endif
      cluIndex += startIndex[entry->non_local_process()];

      // indicate - there is no reason now why a vertex may not
      // match non-locally with a vertex that was matched non-locally
      // with a vertex from another processor during a high-to-low
      // communication for example...otherwise the cluster weight
      // information is redundant?

      for (i = 0; i < numLocals; ++i)
        match_vector_[locals[i] - minimum_vertex_index_] = cluIndex;
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

int model_coarsener_2d::accept(int locVertex, int nonLocCluWt, int highToLow,
                                 int maxWt) {
  int locVertexIndex = locVertex - minimum_vertex_index_;
  int matchValue = match_vector_[locVertexIndex];

  if (matchValue < NON_LOCAL_MATCH) {
    if (cluster_weights_[matchValue] + nonLocCluWt >= maxWt)
      return 0;
    else {
      // ###
      // the non-local requesting vertices will have the same
      // cluster index as locVertex
      // ###

      cluster_weights_[matchValue] += nonLocCluWt;
      return 1;
    }
  } else {

    int nonLocReq = matchValue - NON_LOCAL_MATCH;
    int proc = nonLocReq / (number_of_vertices_ / processors_);

    if ((highToLow && proc < rank_) || (!highToLow && proc > rank_))
      return 0;

    // ###
    // here think about allowing more than one cross-processor match
    // ###

    if (table_->cluster_index(nonLocReq) != -1)
      return 0;

    int cluWt = vertex_weights_[locVertexIndex] + nonLocCluWt;

    if (cluWt >= maxWt)
      return 0;

    else {
      // ###
      // the non-local requesting vertices will have the
      // same cluster index as locVertex
      // remove locVertex from the request table
      // ###

      table_->remove_local(nonLocReq, locVertex, vertex_weights_[locVertexIndex]);
      match_vector_[locVertexIndex] = cluster_index_;
      cluster_weights_[cluster_index_++] = cluWt;

      return 1;
    }
  }
}

void model_coarsener_2d::permute_vertices_arrays(dynamic_array<int> &verts,
                                                 int nLocVerts) {
  int i;
  switch (vertex_visit_order_) {
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
    verts.random_permutation();
    break;

  case INCREASING_WEIGHT_ORDER:
    for (i = 0; i < nLocVerts; ++i) {
      verts[i] = i;
    }
    verts.sort_between_using_another_array(
        0, nLocVerts - 1, vertex_weights_,
        parkway::utility::sort_order::INCREASING);
    break;

  case DECREASING_WEIGHT_ORDER:
    for (i = 0; i < nLocVerts; ++i) {
      verts[i] = i;
    }
    verts.sort_between_using_another_array(
        0, nLocVerts - 1, vertex_weights_,
        parkway::utility::sort_order::DECREASING);
    break;

  default:
    for (i = 0; i < nLocVerts; ++i) {
      verts[i] = i;
    }
    verts.random_permutation();
    break;
  }
}



void model_coarsener_2d::print_visit_order(int variable) const {
  switch (variable) {
  case INCREASING_ORDER:
    info("inc-idx");
    break;

  case DECREASING_ORDER:
    info("dec-idx");
    break;

  case INCREASING_WEIGHT_ORDER:
    info("inc_wt");
    break;

  case DECREASING_WEIGHT_ORDER:
    info("dec-wt");
    break;

  default:
    info("rand");
    break;
  }
}

}  // namespace parallel
}  // namespace parkway
