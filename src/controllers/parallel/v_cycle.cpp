#include "controllers/parallel/v_cycle.hpp"
#include "data_structures/complete_binary_tree.hpp"
#include "utility/logging.hpp"

namespace parkway {
namespace parallel {

v_cycle::v_cycle(
    parallel::restrictive_coarsening &rc, parallel::coarsener &c,
    parallel::refiner &r, parkway::serial::controller &ref, int rank, int nP,
    int percentile, int inc, int approxRef, int limit, double limitAsPercent)
    : parallel::controller(c, r, ref, rank, nP, percentile, inc, approxRef),
      restrictive_coarsening_(rc) {
  if (limit)
    limit_on_cycles_ = limit;
  else
    limit_on_cycles_ = LARGE_CONSTANT;

  if (limitAsPercent)
    limit_as_percent_of_cut_ = limitAsPercent;
  else
    limitAsPercent = 0;

  minimum_inter_vertex_index_ = 0;

  map_to_inter_vertices_.reserve(0);
  map_to_orig_vertices_.reserve(0);
}

v_cycle::~v_cycle() {}

void v_cycle::display_options() const {
  info("|--- PARA_CONTR (# parts = %i)\n"
       "|- VCYCLE:  pRuns = %i kT = %.2f rkT = %.2f appRef = %i wTF = %i",
       total_number_of_parts_, number_of_runs_, keep_partitions_within_,
       reduction_in_keep_threshold_, approximate_refine_,
       write_partition_to_file_);
  print_type();
  info(" lim = %i %%min = %.2f start %%le = %i %%le inc = %i\n|\n",
       limit_on_cycles_, limit_as_percent_of_cut_, start_percentile_,
       percentile_increment_);
}

void v_cycle::set_weight_constraints(MPI_Comm comm) {
#ifdef DEBUG_CONTROLLER
  assert(hgraph);
  assert(numTotalParts != -1);
  assert(balConstraint > 0.0 && balConstraint < 1.0);
#endif

  int locGraphWt;
  int totGraphWt;
  int maxVertWt;

  double avePartWt;

  locGraphWt = hypergraph_->vertex_weight();

  MPI_Allreduce(&locGraphWt, &totGraphWt, 1, MPI_INT, MPI_SUM, comm);

  avePartWt = static_cast<double>(totGraphWt) / total_number_of_parts_;
  maximum_part_weight_ = static_cast<int>(floor(avePartWt + avePartWt *
                                                            balance_constraint_));
  maxVertWt = static_cast<int>(floor(avePartWt * balance_constraint_));

  coarsener_.set_maximum_vertex_weight(maxVertWt);
  restrictive_coarsening_.set_maximum_vertex_weight(maxVertWt);

  coarsener_.set_total_hypergraph_weight(totGraphWt);
  restrictive_coarsening_.set_total_graph_weight(totGraphWt);

  serial_controller_.set_maximum_vertex_weight(maxVertWt);
}

void v_cycle::record_v_cycle_partition(
    const parallel::hypergraph &h, int numIteration) {
#ifdef DEBUG_CONTROLLER
  assert(numIteration >= 0);
#endif

  int i;

  int minVertexIndex;
  int numLocalVertices;

  ds::dynamic_array<int> *bestPartVector;

  auto pVector = h.partition_vector();
  minVertexIndex = h.minimum_vertex_index();
  numLocalVertices = h.number_of_vertices();

  if (numIteration == 0) {
    bestPartVector = new ds::dynamic_array<int>(numLocalVertices);

    for (i = 0; i < numLocalVertices; ++i)
      bestPartVector->at(i) = pVector[i];

    minimum_local_current_vertices_.push(minVertexIndex);
    number_of_local_current_vertices_.push(numLocalVertices);
    best_v_cycle_partition_.push(bestPartVector);
  } else {
    bestPartVector = best_v_cycle_partition_.top();

#ifdef DEBUG_CONTROLLER
    assert(bestPartVector);
#endif

    bestPartVector->reserve(numLocalVertices);
    for (i = 0; i < numLocalVertices; ++i)
      bestPartVector->at(i) = pVector[i];

    minimum_local_current_vertices_.pop();
    number_of_local_current_vertices_.pop();
    minimum_local_current_vertices_.push(minVertexIndex);
    number_of_local_current_vertices_.push(numLocalVertices);
  }

#ifdef DEBUG_CONTROLLER
  assert(minLocCurrVertId.getNumElem() > 0);
  assert(numLocCurrVerts.getNumElem() > 0);
  assert(h.getNumPartitions() == 1);
  assert(bestVCyclePartition.getNumElem() > 0);
#endif
}

void v_cycle::gather_in_v_cycle_partition(parallel::hypergraph &h,
                                                              int cut,
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
  int origV;
  int arrayLen;
  int totToSend;
  int totToRecv;

#ifdef DEBUG_CONTROLLER
  int numTotVertices = h.getNumTotalVertices();
#endif

  auto myPartVector = h.partition_vector();
  int numLocalVertices = h.number_of_vertices();
  auto myVtoOrigV = h.to_origin_vertex();

  minStoredVertexIndex = minimum_local_current_vertices_.top();
  minimum_local_current_vertices_.pop();
  numStoredLocVertices = number_of_local_current_vertices_.top();
  number_of_local_current_vertices_.pop();
  maxStoredVertexIndex = minStoredVertexIndex + numStoredLocVertices;

  ds::dynamic_array<int> *bestPartVector = best_v_cycle_partition_.top();
  best_v_cycle_partition_.pop();

#ifdef DEBUG_CONTROLLER
  assert(bestPartVector);
#endif

  ds::dynamic_array<int> minStoredIDs(processors_);
  ds::dynamic_array<int> localSendArrayVerts;
  ds::dynamic_array<int> procs(processors_);

  MPI_Allgather(&minStoredVertexIndex, 1, MPI_INT, minStoredIDs.data(), 1,
                MPI_INT, comm);

  for (i = 0; i < processors_; ++i)
    procs[i] = i;

  ds::complete_binary_tree<int> storedVtoProc(procs.data(),
                                               minStoredIDs.data(), processors_);

  for (i = 0; i < processors_; ++i)
    send_lens_[i] = 0;

  for (i = 0; i < numLocalVertices; ++i) {
    origV = myVtoOrigV[i];
#ifdef DEBUG_CONTROLLER
    assert(origV >= 0 && origV < numTotVertices);
#endif

    if (origV >= minStoredVertexIndex && origV < maxStoredVertexIndex) {
      myPartVector[i] = bestPartVector->at(origV - minStoredVertexIndex);
    } else {
      j = storedVtoProc.root_value(origV);

      data_out_sets_[j][send_lens_[j]++] = origV;
      data_out_sets_[j][send_lens_[j]++] = i;
    }
  }

  ij = 0;

  for (i = 0; i < processors_; ++i) {
    send_lens_[i] = send_lens_[i] >> 1;
    send_displs_[i] = ij;
    ij += send_lens_[i];
  }

  send_array_.reserve(ij);
  localSendArrayVerts.reserve(ij);
  totToSend = ij;

  j = 0;

  for (i = 0; i < processors_; ++i) {
    arrayLen = send_lens_[i] << 1;
    for (ij = 0; ij < arrayLen;) {
      send_array_[j] = data_out_sets_[i][ij++];
      localSendArrayVerts[j++] = data_out_sets_[i][ij++];
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
    send_array_[i] = bestPartVector->at(receive_array_[i] - minStoredVertexIndex);
  }

  receive_array_.reserve(totToSend);

  MPI_Alltoallv(send_array_.data(), receive_lens_.data(),
                receive_displs_.data(), MPI_INT, receive_array_.data(),
                send_lens_.data(), send_displs_.data(), MPI_INT, comm);

  for (i = 0; i < totToSend; ++i) {
    myPartVector[localSendArrayVerts[i]] = receive_array_[i];
  }

  h.set_cut(0, cut);
}

void v_cycle::project_v_cycle_partition(
    parallel::hypergraph &cG, parallel::hypergraph &fG, MPI_Comm comm) {
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

  int numLocalFineVertices = fG.number_of_vertices();

  if (map_to_inter_vertices_.capacity() != 0) {
      fG.set_number_of_partitions(1);

    int numLocalCoarseVertices = cG.number_of_vertices();
    int numTotalCoarseVertices = cG.total_number_of_vertices();
    int coarseCut = cG.cut(0);

    fG.set_cut(0, coarseCut);

    auto coarsePartVector = cG.partition_vector();
    auto finePartVector = fG.partition_vector();
    auto fineMatVector = fG.match_vector();

    int minCoarseVertexId = cG.minimum_vertex_index();
    int numRequestingLocalVerts;
    int cVertPerProc = numTotalCoarseVertices / processors_;
    int cVertex;
    int vPart;
    int totToRecv;
    int totToSend;
    int sendLength;

    int j;
    int ij;

    ds::dynamic_array<int> interGraphPVector(numLocalCoarseVertices);
    ds::dynamic_array<int> requestingLocalVerts;

    for (i = 0; i < processors_; ++i)
      send_lens_[i] = 0;

    for (i = 0; i < numLocalCoarseVertices; ++i) {
      cVertex = map_to_inter_vertices_[i];
      vPart = coarsePartVector[i];

#ifdef DEBUG_CONTROLLER
      assert(cVertex >= 0 && cVertex < numTotalCoarseVertices);
      assert(vPart >= 0 && vPart < numTotalParts);
#endif

      ij = std::min(cVertex / cVertPerProc, processors_ - 1);

      if (ij == rank_) {
        interGraphPVector[cVertex - minCoarseVertexId] = vPart;
      } else {
        data_out_sets_[ij][send_lens_[ij]++] = cVertex;
        data_out_sets_[ij][send_lens_[ij]++] = vPart;
      }
    }

    ij = 0;

    for (i = 0; i < processors_; ++i) {
      send_displs_[i] = ij;
      ij += send_lens_[i];
    }

    send_array_.resize(ij);
    totToSend = ij;
    ij = 0;

    for (i = 0; i < processors_; ++i) {
      j = 0;
      sendLength = send_lens_[i];
      while (j < sendLength) {
        send_array_[ij++] = data_out_sets_[i][j++];
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
      ij = std::min(cVertex / cVertPerProc, processors_ - 1);

      if (ij == rank_) {
        finePartVector[i] = interGraphPVector[cVertex - minCoarseVertexId];
      } else {
        data_out_sets_[ij][send_lens_[ij]++] = i;
        data_out_sets_[ij][send_lens_[ij]++] = cVertex;
      }
    }

    // ###
    // prepare to send requests for partition vector values
    // ###

    ij = 0;

    for (i = 0; i < processors_; ++i) {
      send_displs_[i] = ij;
      ij += send_lens_[i] >> 1;
    }

    send_array_.reserve(ij);
    requestingLocalVerts.reserve(ij);

    numRequestingLocalVerts = ij;
    totToSend = ij;
    ij = 0;

    for (i = 0; i < processors_; ++i) {
      j = 0;
      sendLength = send_lens_[i];
      while (j < sendLength) {
        requestingLocalVerts[ij] = data_out_sets_[i][j++];
        send_array_[ij++] = data_out_sets_[i][j++];
      }

      send_lens_[i] = sendLength >> 1;
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
      fG.project_partitions(cG, comm);
  }

  map_to_inter_vertices_.reserve(numLocalFineVertices);
  minimum_inter_vertex_index_ = fG.minimum_vertex_index();

  for (i = 0; i < numLocalFineVertices; ++i)
    map_to_inter_vertices_[i] = minimum_inter_vertex_index_ + i;
}

void v_cycle::shuffle_v_cycle_vertices_by_partition(
    parallel::hypergraph &h, MPI_Comm comm) {
  // function spec:
  // as the vertices are shuffled modify the mapToInterVerts
  // check that the sum of the lengths of the mapToInterVerts
  // across the processors is equal to the number of totalVertices
  // of h
  //

  int i;
  int j;
  int ij;

  int minLocVertIdBefShuff = h.minimum_vertex_index();

#ifdef DEBUG_CONTROLLER
  int numLocVertBefShuff = h.getNumLocalVertices();
  int maxLocVertIdBefShuff = minLocVertIdBefShuff + numLocVertBefShuff;
#endif

  int numTotVertices = h.total_number_of_vertices();
  int numVBefPerProc = numTotVertices / processors_;

#ifdef DEBUG_CONTROLLER
  assert(minInterVertIndex == minLocVertIdBefShuff);
  assert(numLocVertBefShuff == mapToInterVerts.capacity());
#endif

    h.shuffle_vertices_by_partition(total_number_of_parts_, comm);

  int numLocVertAftShuff = h.number_of_vertices();
  int minLocVertIdAftShuff = h.minimum_vertex_index();
  int totToSend;
  int totToRecv;
  int sendLength;
  int numRequestingLocalVerts;

  auto vToOrigV = h.to_origin_vertex();

  dynamic_array<int> newToInter(numLocVertAftShuff);
  dynamic_array<int> requestingLocalVerts;

  for (i = 0; i < processors_; ++i)
    send_lens_[i] = 0;

  for (i = 0; i < numLocVertAftShuff; ++i) {
    j = vToOrigV[i];
    ij = std::min(j / numVBefPerProc, processors_ - 1);

    if (ij == rank_) {
      newToInter[i] = map_to_inter_vertices_[j - minimum_inter_vertex_index_];
    } else {
      data_out_sets_[ij][send_lens_[ij]++] = i;
      data_out_sets_[ij][send_lens_[ij]++] = j;
    }
  }

  ij = 0;
  for (i = 0; i < processors_; ++i) {
    send_displs_[i] = ij;
    ij += send_lens_[i] >> 1;
  }

  send_array_.resize(ij);
  requestingLocalVerts.resize(ij);
  numRequestingLocalVerts = ij;
  totToSend = ij;
  ij = 0;

  for (i = 0; i < processors_; ++i) {
    j = 0;
    sendLength = send_lens_[i];
    while (j < sendLength) {
      requestingLocalVerts[ij] = data_out_sets_[i][j++];
      send_array_[ij++] = data_out_sets_[i][j++];
    }

    send_lens_[i] = send_lens_[i] >> 1;
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

  receive_array_.resize(ij);
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
    send_array_[i] = map_to_inter_vertices_[j];
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

  for (i = 0; i < numRequestingLocalVerts; ++i) {
    j = requestingLocalVerts[i];
#ifdef DEBUG_CONTROLLER
    assert(j >= 0 && j < numLocVertAftShuff);
    assert(receive_array_[i] >= 0 && receive_array_[i] < numTotVertices);
#endif
    newToInter[j] = receive_array_[i];
  }

  minimum_inter_vertex_index_ = minLocVertIdAftShuff;
  map_to_inter_vertices_ = newToInter;
}

void v_cycle::shift_v_cycle_vertices_to_balance(
    parallel::hypergraph &h, MPI_Comm comm) {
  // function spec:
  // as the vertices are shuffled modify the mapToInterVerts
  // check that the sum of the lengths of the mapToInterVerts
  // across the processors is equal to the number of totalVertices
  // of h
  //

  int i;
  int j;

  int numLocalVertices = h.number_of_vertices();
  int minVertexIndex = h.minimum_vertex_index();
  int maxVertexIndex = minVertexIndex + numLocalVertices;
  int numTotVertices = h.total_number_of_vertices();
  int vPerProc = numTotVertices / processors_;
  int numMyNewVertices;

  if (rank_ != processors_ - 1)
    numMyNewVertices = vPerProc;
  else
    numMyNewVertices = vPerProc + (numTotVertices % processors_);

  dynamic_array<int> minNewIndex(processors_);
  dynamic_array<int> maxNewIndex(processors_);

#ifdef DEBUG_CONTROLLER
  assert(minVertexIndex == minInterVertIndex);
#endif

    h.shift_vertices_to_balance(comm);

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
        std::max(numLocalVertices - (std::max(maxVertexIndex - maxNewIndex[i], 0) +
                                std::max(minNewIndex[i] - minVertexIndex, 0)),
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

  MPI_Alltoallv(map_to_inter_vertices_.data(), send_lens_.data(),
                send_displs_.data(), MPI_INT, newToInter->data(),
                receive_lens_.data(), receive_displs_.data(), MPI_INT, comm);

  minimum_inter_vertex_index_ = minNewIndex[rank_];
  map_to_inter_vertices_ = *newToInter;
}

void v_cycle::update_map_to_orig_vertices(MPI_Comm comm) {
#ifdef DEBUG_CONTROLLER
  int numLocalVertices = hgraph->getNumLocalVertices();
#endif

  int minLocVertIndex = hypergraph_->minimum_vertex_index();
  int numTotalVertices = hypergraph_->total_number_of_vertices();
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

  dynamic_array<int> newMapToOrig(number_of_orig_local_vertices_);
  dynamic_array<int> copyOfSendArray;

  for (i = 0; i < processors_; ++i)
    send_lens_[i] = 0;

  for (i = 0; i < number_of_orig_local_vertices_; i++) {
    vertex = map_to_inter_vertices_[i];

#ifdef DEBUG_CONTROLLER
    assert(vertex >= 0 && vertex < numTotalVertices);
#endif

    ij = std::min(vertex / vertPerProc, processors_ - 1);

    if (ij == rank_) {
      newMapToOrig[i] = map_to_orig_vertices_[vertex - minLocVertIndex];
    } else {
      data_out_sets_[ij][send_lens_[ij]++] = i;
      data_out_sets_[ij][send_lens_[ij]++] = vertex;
    }
  }

  ij = 0;

  for (i = 0; i < processors_; ++i) {
    send_displs_[i] = ij;
    ij += send_lens_[i] >> 1;
  }

  send_array_.reserve(ij);
  copyOfSendArray.reserve(ij);
  totToSend = ij;
  ij = 0;

  for (i = 0; i < processors_; ++i) {
    j = 0;
    sendLength = send_lens_[i];
    while (j < sendLength) {
      copyOfSendArray[ij] = data_out_sets_[i][j++];
      send_array_[ij++] = data_out_sets_[i][j++];
    }

    send_lens_[i] = send_lens_[i] >> 1;
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
    send_array_[i] = map_to_orig_vertices_[receive_array_[i] - minLocVertIndex];
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
    newMapToOrig[copyOfSendArray[i]] = receive_array_[i];
  }

#ifdef DEBUG_CONTROLLER
  for (i = 0; i < numOrigLocVerts; ++i)
    assert(newMapToOrig[i] >= 0 && newMapToOrig[i] < numTotalVertices);
#endif

  map_to_orig_vertices_ = newMapToOrig;
}

void v_cycle::reset_structures() {
    hypergraph_->reset_vectors();

  map_to_inter_vertices_.reserve(0);

#ifdef DEBUG_CONTROLLER
  assert(hgraphs.getNumElem() == 0);
  assert(bestVCyclePartition.getNumElem() == 0);
#endif
}

}  // namespace parallel
}  // namespace parkway
