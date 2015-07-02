// ### ParaApproxFCCoarsener.cpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 31/12/2004: Last Modified
//
// ###

#include "coarseners/parallel/approximate_first_choice_coarsener.hpp"
#include "data_structures/map_to_pos_int.hpp"
#include "data_structures/internal/table_utils.hpp"
#include "utility/logging.hpp"
#include "Macros.h"
#include <iostream>

namespace parkway {
namespace parallel {

approximate_first_choice_coarsener::approximate_first_choice_coarsener(
    int rank, int nProcs, int nParts, int percentile, int inc, int vertVisOrder,
    int matchReqOrder, int divByWt, int divByLen)
    : approximate_coarsener(rank, nProcs, nParts, percentile, inc) {
  vertexVisitOrder = vertVisOrder;
  matchRequestVisitOrder = matchReqOrder;
  divByCluWt = divByWt;
  divByHedgeLen = divByLen;
  limitOnIndexDuringCoarsening = 0;

  table = nullptr;
}

void approximate_first_choice_coarsener::dispCoarseningOptions() const {
  info("|--- PARA_C:\n"
       "|- ApproxPFC: r = %.2f min = %i vvo = ",
       reduction_ratio_, minimum_number_of_nodes_);
  printVisitOrder(vertexVisitOrder);
  info(" mvo = ");
  printVisitOrder(matchRequestVisitOrder);
  info(" divWt = %i divLen = %i percentile = %i inc = %i\n|\n",
       divByCluWt, divByHedgeLen, startPercentile, increment);

}

void approximate_first_choice_coarsener::buildAuxiliaryStructs(
    int numTotPins, double aveVertDeg, double aveHedgeSize) {
  // ###
  // build the ds::match_request_table
  // ###

  int i = static_cast<int>(ceil(static_cast<double>(numTotPins) / aveVertDeg));

  // table = new ds::match_request_table(Funct::setTableSize(i/processors_));
  table = new ds::match_request_table(ds::internal::table_utils::table_size(i / processors_));
  // i = Shiftl(static_cast<int>(ceil(aveVertDeg*aveHedgeSize)),4);

  // ###
  // build the ConnVertTable
  // ###

  // connTable = new ConnVertTable(Funct::setTableSize(i));
}

void approximate_first_choice_coarsener::release_memory() {
  hyperedge_weights_.reserve(0);
  hyperedge_offsets_.reserve(0);
  local_pin_list_.reserve(0);

  vertex_to_hyperedges_offset_.reserve(0);
  vertex_to_hyperedges_.reserve(0);
  allocated_hyperedges_.reserve(0);
  cluster_weights_.reserve(0);

    free_memory();
}

parallel::hypergraph *approximate_first_choice_coarsener::coarsen(
    parallel::hypergraph &h, MPI_Comm comm) {
  load(h, comm);

  if (number_of_vertices_ < minimum_number_of_nodes_ || h.dont_coarsen()) {
    currPercentile = startPercentile;
    return nullptr;
  }

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
  // ConnVertData *vData;
  // ConnVertData **vDataArray;

  dynamic_array<int> neighVerts;
  dynamic_array<int> neighCluWts;
  dynamic_array<double> connectVals;
  dynamic_array<int> vertices(number_of_local_vertices_);

  permuteVerticesArray(vertices, number_of_local_vertices_);

  if (parkway::utility::status::handler::progress_enabled()) {
    for (i = 0; i < number_of_local_vertices_; ++i) {
      if (vertex_weights_[i] > maxLocWt)
        maxLocWt = vertex_weights_[i];
    }

    MPI_Reduce(&maxLocWt, &maxWt, 1, MPI_INT, MPI_MAX, 0, comm);

    info("[PFCC] %i %i %i ", maximum_vertex_weight_, maxWt, aveVertexWt);
  }

  metric = static_cast<double>(number_of_local_vertices_) / reduction_ratio_;
  limitOnIndexDuringCoarsening =
      number_of_local_vertices_ - static_cast<int>(floor(metric - 1.0));
  cluster_index_ = 0;
  stop_coarsening_ = 0;

  for (; index < number_of_local_vertices_; ++index) {
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
            // ###
            // now try not to check the weight before checking the vertex
            // ###

            auto neighbourLoc = matchInfoLoc.get_careful(candidatV);
            // vData = connTable->getDataStruct(candidatV);

            if (neighbourLoc >= 0) {
              if (divByHedgeLen)
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
                              table->cluster_weight(nonLocV) + aveVertexWt;
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

                candVwt = table->cluster_weight(candidatV);

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

              if (divByHedgeLen)
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

        if (pairWt <= maximum_vertex_weight_) {
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

           if(divide_by_cluster_weight_)
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

      matchInfoLoc.clear();
      numVisited = 0;
      // connTable->clearTable();

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
              table->add_local(nonLocV, vertex + minimum_vertex_index_,
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

          table->add_local(bestMatch, vertex + minimum_vertex_index_, vertex_weights_[vertex],
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

      if (index > limitOnIndexDuringCoarsening) {
        reducedBy = static_cast<double>(number_of_local_vertices_) /
                    (numNotMatched + cluster_index_ + table->size());
        break;
      }
    }
  }

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
    setRequestArrays(i);
      send_from_data_out(comm); // actually sending requests
    setReplyArrays(i, maximum_vertex_weight_);
      send_from_data_out(comm); // actually sending replies
    processReqReplies();
  }

  set_cluster_indices(comm);

  if (static_cast<double>(number_of_vertices_) / total_number_of_clusters_ <
      MIN_ALLOWED_REDUCTION_RATIO) {
    stop_coarsening_ = 1;
  }

  table->clear();

  // ###
  // now construct the coarse hypergraph using the matching vector
  // ###

  return (contract_hyperedges(h, comm));
}

void approximate_first_choice_coarsener::setRequestArrays(int highToLow) {
  int numRequests = table->size();
  int nonLocVertex;
  int cluWt;

  int i;
  int procRank;

  ds::match_request_table::entry *entry;
  ds::dynamic_array<ds::match_request_table::entry *> entryArray = table->get_entries();

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

    if ((cluWt > 0) && ((highToLow && procRank < rank_) ||
                        (!highToLow && procRank > rank_))) {
      data_out_sets_[procRank][send_lens_[procRank]++] = nonLocVertex;
      data_out_sets_[procRank][send_lens_[procRank]++] = cluWt;
    }
  }
}

void approximate_first_choice_coarsener::setReplyArrays(int highToLow, int maxVWt) {
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

    if (matchRequestVisitOrder == RANDOM_ORDER) {
      visitOrderLen = receive_lens_[i] >> 1;
      visitOrder.reserve(visitOrderLen);

      for (j = 0; j < visitOrderLen; ++j)
        visitOrder[j] = j;

      // ###
      // processing match requests in random order
      // ###

      Funct::randomPermutation(visitOrder.data(), visitOrderLen);

      for (l = 0; l < visitOrderLen; ++l) {

        j = visitOrder[l] << 1;

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

void approximate_first_choice_coarsener::processReqReplies() {
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
        entry = table->get_entry(vNonLocReq);

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

        entry = table->get_entry(vNonLocReq);
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

void approximate_first_choice_coarsener::permuteVerticesArray(dynamic_array<int> &verts, int nLocVerts) {
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

void approximate_first_choice_coarsener::setClusterIndices(MPI_Comm comm) {
  dynamic_array<int> numClusters(processors_);
  dynamic_array<int> startIndex(processors_);

  MPI_Allgather(&cluster_index_, 1, MPI_INT, numClusters.data(), 1, MPI_INT,
                comm);

  int index = 0;
  int i;

  ds::match_request_table::entry *entry;

  int numLocals;
  int cluIndex;
  int numEntries = table->size();

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

  ds::dynamic_array<ds::match_request_table::entry *> entryArray = table->get_entries();

  for (index = 0; index < numEntries; ++index) {
    entry = entryArray[index];

#ifdef DEBUG_COARSENER
    assert(entry);
#endif

    cluIndex = entry->cluster_index();

    if (cluIndex >= 0 && cluIndex != MATCHED_LOCALLY) {
      numLocals = entry->number_local();
      auto locals = entry->local_vertices_array();

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
    assert(matchVector[index] >= 0 && matchVector[index] < totalVertices);
#endif
}

int approximate_first_choice_coarsener::accept(int locVertex, int nonLocCluWt, int highToLow,
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

    if (table->cluster_index(nonLocReq) != -1)
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

      table->remove_local(nonLocReq, locVertex, vertex_weights_[locVertexIndex]);
      match_vector_[locVertexIndex] = cluster_index_;
      cluster_weights_[cluster_index_++] = cluWt;

      return 1;
    }
  }
}

void approximate_first_choice_coarsener::printVisitOrder(int variable) const {
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
