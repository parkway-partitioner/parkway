// ### FMRefiner.cpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 30/11/2004: Last Modified
//
// ###

#include "refiners/serial/fm_refiner.hpp"

namespace parkway {
namespace serial {

fm_refiner::fm_refiner(int max, int insMethod, int ee, int dL) : refiner(dL) {
  bucket_arrays_length_ = 0;
  maximum_possible_gain_ = 0;
  maximum_non_positive_moves_ = 0;
  buckets_.reserve(0);
  move_list_.reserve(0);
  vertex_gains_.reserve(0);
  vertex_in_part_.reserve(0);
  loaded_.reserve(0);

  ee_threshold_ = ee;
  insert_method_ = insMethod;
  maximum_part_weight_ = max;
  number_of_parts_ = 2;
  part_weights_.reserve(2);
  bucket_arrays_.reserve(2);
  number_of_buckets_in_array_.reserve(2);
  maximum_gain_entries_.reserve(2);
  maximum_gains_.reserve(2);

  bucket_arrays_[0] = nullptr;
  bucket_arrays_[1] = nullptr;
}

fm_refiner::~fm_refiner() {
  dynamic_memory::delete_pointer<NodeArray>(bucket_arrays_[0]);
  dynamic_memory::delete_pointer<NodeArray>(bucket_arrays_[1]);
}

void fm_refiner::display_options(std::ostream &out) const {
  switch (dispOption) {
  case SILENT:
    break;

  default:

    out << "|- FM:"
        << " qDis = ";
    printQdis(out);
    out << " eeT = " << ee_threshold_ << " acc = " << accept_proportion_ << std::endl
        << "|" << std::endl;
    break;
  }
}

void fm_refiner::printQdis(std::ostream &out) const {
  switch (insert_method_) {
  case FIFO:
    out << "FIFO";
    break;

  default:
    out << "LIFO";
    break;
  }
}

void fm_refiner::build_buckets() {
  int i;
  int j;

  int maxVertDeg = 0;
  int maxHedgeWt = 0;

  buckets_.reserve(numVertices);
  vertex_gains_.reserve(numVertices);
  loaded_.reserve(numVertices);
  move_list_.reserve(numVertices);
  vertex_in_part_.reserve(Shiftl(numHedges, 1));

  for (i = 0; i < numVertices; ++i) {
    buckets_[i] = new bucket_node;
    buckets_[i]->vertex_id = i;
    buckets_[i]->previous = nullptr;
    buckets_[i]->next = nullptr;

    j = vOffsets[i + 1] - vOffsets[i];

    if (j > maxVertDeg)
      maxVertDeg = j;
  }

  for (i = 0; i < numHedges; ++i) {
    if (hEdgeWeight[i] > maxHedgeWt) {
      maxHedgeWt = hEdgeWeight[i];
    }
  }

  maximum_possible_gain_ = maxVertDeg * maxHedgeWt;
  bucket_arrays_length_ = Or(Shiftl(maximum_possible_gain_, 1), 0x1);

  bucket_arrays_[0] = new NodeArray(bucket_arrays_length_);
  bucket_arrays_[1] = new NodeArray(bucket_arrays_length_);

  for (i = 0; i < bucket_arrays_length_; ++i) {
    (*bucket_arrays_[0])[i].next = nullptr;
    (*bucket_arrays_[0])[i].vertex_id = -1;
    (*bucket_arrays_[1])[i].next = nullptr;
    (*bucket_arrays_[1])[i].vertex_id = -1;
  }
}

void fm_refiner::restore_buckets() {
  int i;

  for (i = 0; i < bucket_arrays_length_; ++i) {
    (*bucket_arrays_[0])[i].next = nullptr;
    (*bucket_arrays_[0])[i].vertex_id = -1;
    (*bucket_arrays_[1])[i].next = nullptr;
    (*bucket_arrays_[1])[i].vertex_id = -1;
  }

  for (i = 0; i < numVertices; ++i) {
    buckets_[i]->previous = nullptr;
    buckets_[i]->next = nullptr;
  }
}

void fm_refiner::destroy_buckets() {
  int i;

  for (i = 0; i < numVertices; ++i)
    dynamic_memory::delete_pointer<bucket_node>(buckets_[i]);

  for (i = 0; i < 2; ++i)
    dynamic_memory::delete_pointer<NodeArray>(bucket_arrays_[i]);
}

void fm_refiner::initialise_partition_structure() {
  int i;
  int j;
  int hEdgeOffset;

  int endOffset;

  part_weights_[0] = 0;
  part_weights_[1] = 0;

  for (i = 0; i < numVertices; ++i) {
    part_weights_[partition_vector_[i]] += vWeight[i];
  }

  for (i = 0; i < numHedges; ++i) {
    hEdgeOffset = Shiftl(i, 1);

    vertex_in_part_[hEdgeOffset] = 0;
    vertex_in_part_[Or(hEdgeOffset, 0x1)] = 0;

    endOffset = hEdgeOffsets[i + 1];

    for (j = hEdgeOffsets[i]; j < endOffset; ++j) {
      ++vertex_in_part_[Or(hEdgeOffset, partition_vector_[pinList[j]])];
    }
  }
}

void fm_refiner::prepare_vertex_gains() {
  int i;
  int j;
  int endOffset;
  int hEdgeOffset;

  int vPart;
  int hEdge;

  // ###
  // compute the vertex gains
  // ###

  for (i = 0; i < 2; ++i) {
    maximum_gains_[i] = -maximum_possible_gain_;
    maximum_gain_entries_[i] = 0;
    number_of_buckets_in_array_[i] = 0;
  }

  for (i = 0; i < numVertices; ++i) {
    endOffset = vOffsets[i + 1];
    vPart = partition_vector_[i];

#ifdef DEBUG_FM_REFINER
    assert(vPart == 0 || vPart == 1);
#endif

    vertex_gains_[i] = 0;

    for (j = vOffsets[i]; j < endOffset; ++j) {
      hEdge = vToHedges[j];
      hEdgeOffset = Shiftl(hEdge, 1);

      if (vertex_in_part_[Or(hEdgeOffset, vPart)] == 1)
        vertex_gains_[i] += hEdgeWeight[hEdge];

      if (vertex_in_part_[Or(hEdgeOffset, Xor(vPart, 1))] == 0)
        vertex_gains_[i] -= hEdgeWeight[hEdge];
    }

    move_to_bucket_array(vPart, vertex_gains_[i], i);
    ++number_of_buckets_in_array_[vPart];
  }
}

void fm_refiner::initialise_gains_1_to_0() {
  int i;
  int j;

  int hEdge;
  int endOffset;

  // ###
  // compute vertex gains in 1 to 0 direction
  // ###

  for (i = 0; i < 2; ++i) {
    maximum_gains_[i] = -maximum_possible_gain_;
    maximum_gain_entries_[i] = 0;
    number_of_buckets_in_array_[i] = 0;
  }

  for (i = 0; i < numVertices; ++i) {
    if (partition_vector_[i] == 1) {
      endOffset = vOffsets[i + 1];
      vertex_gains_[i] = 0;

      for (j = vOffsets[i]; j < endOffset; ++j) {
        hEdge = vToHedges[j];

        if (vertex_in_part_[Shiftl(hEdge, 1)] == 0)
          vertex_gains_[i] -= hEdgeWeight[hEdge];

        if (vertex_in_part_[Or(Shiftl(hEdge, 1), 0x1)] == 1)
          vertex_gains_[i] += hEdgeWeight[hEdge];
      }

      move_to_bucket_array(1, vertex_gains_[i], i);

#ifdef DEBUG_FM_REFINER
      assert(buckets[i]->prev != nullptr);
#endif
      ++number_of_buckets_in_array_[1];
    }
  }
}

void fm_refiner::initialise_gains_0_to_1() {
  int i;
  int j;

  int hEdge;
  int endOffset;

  // ###
  // compute vertex gains in 0 to 1 direction
  // ###

  for (i = 0; i < 2; ++i) {
    maximum_gains_[i] = -maximum_possible_gain_;
    maximum_gain_entries_[i] = 0;
    number_of_buckets_in_array_[i] = 0;
  }

  for (i = 0; i < numVertices; ++i) {
    if (partition_vector_[i] == 0) {
      endOffset = vOffsets[i + 1];
      vertex_gains_[i] = 0;

      for (j = vOffsets[i]; j < endOffset; ++j) {
        hEdge = vToHedges[j];

        if (vertex_in_part_[Or(Shiftl(hEdge, 1), 0x1)] == 0)
          vertex_gains_[i] -= hEdgeWeight[hEdge];

        if (vertex_in_part_[Shiftl(hEdge, 1)] == 1)
          vertex_gains_[i] += hEdgeWeight[hEdge];
      }

      move_to_bucket_array(0, vertex_gains_[i], i);

#ifdef DEBUG_FM_REFINER
      assert(buckets[i]->prev != nullptr);
#endif
      ++number_of_buckets_in_array_[0];
    }
  }
}

void fm_refiner::remove_buckets_from_1() {
  int i;
  bucket_node *b;

  bucket_node *array = bucket_arrays_[1]->data();

  int bucketArrayLen = Or(Shiftl(maximum_possible_gain_, 1), 0x1);

  for (i = 0; i < bucketArrayLen; ++i) {
    while (array[i].next != nullptr) {
      b = array[i].next;
#ifdef DEBUG_FM_REFINER
      assert(b);
#endif
      array[i].next = b->next;

      b->previous = nullptr;
      b->next = nullptr;
    }
  }

  number_of_buckets_in_array_[1] = 0;
}

void fm_refiner::remove_buckets_from_0() {
  int i;
  bucket_node *b;

  bucket_node *array = bucket_arrays_[0]->data();

  int bucketArrayLen = Or(Shiftl(maximum_possible_gain_, 1), 0x1);

  for (i = 0; i < bucketArrayLen; ++i) {
    while (array[i].next != nullptr) {
      b = array[i].next;
#ifdef DEBUG_FM_REFINER
      assert(b);
#endif
      array[i].next = b->next;

      b->previous = nullptr;
      b->next = nullptr;
    }
  }

  number_of_buckets_in_array_[0] = 0;
}

void fm_refiner::remove_unlocked_from_bucket_arrays() {
  int i;
  bucket_node *b;

  for (i = 0; i < numVertices; ++i) {
    if (!loaded_(i)) {
      b = buckets_[i];

#ifdef DEBUG_FM_REFINER
      assert(b);
      assert(b->prev != nullptr);
#endif

      b->previous->next = b->next;

      if (b->next)
        b->next->previous = b->previous;

      b->previous = nullptr;
      b->next = nullptr;
    }
  }

  number_of_buckets_in_array_[0] = 0;
  number_of_buckets_in_array_[1] = 0;
}

void fm_refiner::move_to_bucket_array(int vPart, int vGain, int v) {
#ifdef DEBUGH_REFINER
  assert(vPart == 0 || vPart == 1);
  assert(v >= 0 && v < numVertices);
#endif

  int arrayIndex = vGain + maximum_possible_gain_;

  bucket_node *b = buckets_[v];
  bucket_node *bucket;

  // ###
  // first remove from bucket data_
  // ###

  if (b->previous) {
    b->previous->next = b->next;

    if (b->next) {
      b->next->previous = b->previous;
    }
  }

  // ###
  // now insert into appropriate slot
  // ###

  if (vGain > maximum_gains_[vPart]) {
    maximum_gain_entries_[vPart] = arrayIndex;
    maximum_gains_[vPart] = vGain;
  }

  bucket = &(*bucket_arrays_[vPart])[arrayIndex];

  switch (insert_method_) {
  case FIFO: {
    while (bucket->next)
      bucket = bucket->next;

    bucket->next = b;
    b->previous = bucket;
    b->next = nullptr;

    break;
  }
  case LIFO: {
    b->next = bucket->next;
    b->previous = bucket;
    bucket->next = b;

    if (b->next)
      b->next->previous = b;

    break;
  }
  default:
    // ###
    // LIFO
    // ###
    {
      b->next = bucket->next;
      b->previous = bucket;
      bucket->next = b;

      if (b->next)
        b->next->previous = b;

#ifdef DEBUG_FM_REFINER
      assert(bucket->vertex_id == -1);
      assert(bucket->next);
#endif
      break;
    }
  }

  if (maximum_gains_[vPart] > vGain) {
    int index = maximum_gain_entries_[vPart];

    while ((*bucket_arrays_[vPart])[index].next == nullptr) {
#ifdef DEBUG_FM_REFINER
      assert(index >= 0);
#endif
      --index;
    }

    maximum_gain_entries_[vPart] = index;
    maximum_gains_[vPart] = index - maximum_possible_gain_;

#ifdef DEBUG_FM_REFINER
    assert(index >= 0);
    assert((*bucketArrays[vPart])[index].next);
#endif
  }
}

void fm_refiner::remove_from_bucket_array(int v, int vPart, int vGain) {
#ifdef DEBUG_FM_REFINER
  assert(vPart == 0 || vPart == 1);
  assert(v >= 0 && v < numVertices);
#endif

  bucket_node *b = buckets_[v];

#ifdef DEBUG_FM_REFINER
  assert(b->prev);
#endif

  b->previous->next = b->next;

  if (b->next)
    b->next->previous = b->previous;

  b->previous = nullptr;
  b->next = nullptr;

  --number_of_buckets_in_array_[vPart];

#ifdef DEBUG_FM_REFINER
  assert(numBucketsInArray[vPart] >= 0);
#endif

  if (maximum_gains_[vPart] == vGain) {
    int index = vGain + maximum_possible_gain_;

    while ((*bucket_arrays_[vPart])[index].next == nullptr && index > 0)
      --index;

    maximum_gain_entries_[vPart] = index;
    maximum_gains_[vPart] = index - maximum_possible_gain_;

#ifdef DEBUG_FM_REFINER
    assert(index >= 0);
#endif
  }
}

void fm_refiner::adjust_maximum_gain_pointer(int dP) {
  if (number_of_buckets_in_array_[dP] > 0) {
    int index = Shiftl(maximum_possible_gain_, 1);

    while ((*bucket_arrays_[dP])[index].next == nullptr) {
#ifdef DEBUG_FM_REFINER
      assert(index > 0);
#endif
      index--;
    }

    maximum_gain_entries_[dP] = index;
    maximum_gains_[dP] = index - maximum_possible_gain_;

#ifdef DEBUG_FM_REFINER
    assert(index >= 0);
    assert((*bucketArrays[dP])[index].next);
#endif
  }
}

void fm_refiner::update_gains(int v) {
#ifdef DEBUG_FM_REFINER
  assert(v >= 0 && v < numVertices);
#endif

  int i;
  int j;
  int adjVertex;
  int hEdge;
  int endOffset;
  int endHedgeOffset;
  int sourceOffset;
  int destOffset;
  int sP = partition_vector_[v];
  int dP = sP ^ 0x1;

  part_weights_[sP] -= vWeight[v];
  part_weights_[dP] += vWeight[v];

  partition_vector_[v] = dP;

  endOffset = vOffsets[v + 1];

  for (i = vOffsets[v]; i < endOffset; ++i) {
    hEdge = vToHedges[i];
    sourceOffset = Or(Shiftl(hEdge, 1), sP);
    destOffset = Or(Shiftl(hEdge, 1), dP);

    // ###
    // now if the hyperedge has no vertices in
    // destination part, increment the gains of
    // other free vertices in the hyperedge
    // ###

    if (vertex_in_part_[destOffset] == 0) {
      endHedgeOffset = hEdgeOffsets[hEdge + 1];

      for (j = hEdgeOffsets[hEdge]; j < endHedgeOffset; ++j) {
        adjVertex = pinList[j];

#ifdef DEBUG_FM_REFINER
        if (adjVertex != v)
          assert(partitionVector[adjVertex] == sP);
#endif
        if (adjVertex != v && !loaded_(adjVertex)) {
          vertex_gains_[adjVertex] += hEdgeWeight[hEdge];
          move_to_bucket_array(sP, vertex_gains_[adjVertex], adjVertex);
        }
      }
    }

    // ###
    // if the hyperedge has one vertex in the
    // destination part and it is free, decrement
    // its gain
    // ###

    if (vertex_in_part_[destOffset] == 1) {
      endHedgeOffset = hEdgeOffsets[hEdge + 1];

      for (j = hEdgeOffsets[hEdge]; j < endHedgeOffset; ++j) {
        adjVertex = pinList[j];

        if (adjVertex != v && partition_vector_[adjVertex] == dP &&
            !loaded_(adjVertex)) {
          vertex_gains_[adjVertex] -= hEdgeWeight[hEdge];
          move_to_bucket_array(dP, vertex_gains_[adjVertex], adjVertex);
          break;
        }
      }
    }

    // ###
    // move v from sP to dP for hyperedge
    // ###

    --vertex_in_part_[sourceOffset];
    ++vertex_in_part_[destOffset];

    // ###
    // if the hyperedge has no vertices in the sP
    // decrement gains of all free vertices
    // ###

    if (vertex_in_part_[sourceOffset] == 0) {
      endHedgeOffset = hEdgeOffsets[hEdge + 1];

      for (j = hEdgeOffsets[hEdge]; j < endHedgeOffset; ++j) {
        adjVertex = pinList[j];

#ifdef DEBUG_FM_REFINER
        if (adjVertex != v)
          assert(partitionVector[adjVertex] == dP);
#endif

        if (adjVertex != v && !loaded_(adjVertex)) {
          vertex_gains_[adjVertex] -= hEdgeWeight[hEdge];
          move_to_bucket_array(dP, vertex_gains_[adjVertex], adjVertex);
        }
      }
    }

    // ###
    // if the hyperedge has one vertex in the sP
    // increment the gain of this vertex if free
    // ###

    if (vertex_in_part_[sourceOffset] == 1) {
      endHedgeOffset = hEdgeOffsets[hEdge + 1];

      for (j = hEdgeOffsets[hEdge]; j < endHedgeOffset; ++j) {
        adjVertex = pinList[j];

        if (adjVertex != v && partition_vector_[adjVertex] == sP &&
            !loaded_(adjVertex)) {
          vertex_gains_[adjVertex] += hEdgeWeight[hEdge];
          move_to_bucket_array(sP, vertex_gains_[adjVertex], adjVertex);
          break;
        }
      }
    }
  }
}

void fm_refiner::update_gains_1_to_0(int v) {
#ifdef DEBUG_FM_REFINER
  assert(v >= 0 && v < numVertices);
  assert(partitionVector[v] == 1);
#endif

  int i;
  int j;
  int adjVertex;

  int hEdge;
  int endOffset;
  int destOffset;
  int sourceOffset;
  int endHedgeOffset;

  part_weights_[1] -= vWeight[v];
  part_weights_[0] += vWeight[v];

  partition_vector_[v] = 0;

  endOffset = vOffsets[v + 1];

  for (i = vOffsets[v]; i < endOffset; ++i) {
    hEdge = vToHedges[i];
    sourceOffset = Or(Shiftl(hEdge, 1), 0x1);
    destOffset = Shiftl(hEdge, 1);

    // ###
    // now if the hyperedge has no vertices in
    // destination part, increment the gains of
    // other free vertices in the hyperedge
    // ###

    if (vertex_in_part_[destOffset] == 0) {
      endHedgeOffset = hEdgeOffsets[hEdge + 1];

      for (j = hEdgeOffsets[hEdge]; j < endHedgeOffset; ++j) {
        adjVertex = pinList[j];

#ifdef DEBUG_FM_REFINER
        assert(adjVertex >= 0 && adjVertex < numVertices);
#endif
        if (adjVertex != v && partition_vector_[adjVertex] == 1) {
          vertex_gains_[adjVertex] += hEdgeWeight[hEdge];
          move_to_bucket_array(1, vertex_gains_[adjVertex], adjVertex);
        }
      }
    }

    // ###
    // move v from sP to dP for hyperedge
    // ###

    --vertex_in_part_[sourceOffset];
    ++vertex_in_part_[destOffset];

    // ###
    // if the hyperedge has one vertex in the sP
    // increment the gain of this vertex if free
    // ###

    if (vertex_in_part_[sourceOffset] == 1) {
      endHedgeOffset = hEdgeOffsets[hEdge + 1];

      for (j = hEdgeOffsets[hEdge]; j < endHedgeOffset; ++j) {
        adjVertex = pinList[j];

#ifdef DEBUG_REFINER
        assert(adjVertex >= 0 && adjVertex < numVertices);
#endif
        if (adjVertex != v && partition_vector_[adjVertex] == 1) {
          vertex_gains_[adjVertex] += hEdgeWeight[hEdge];
          move_to_bucket_array(1, vertex_gains_[adjVertex], adjVertex);
          break;
        }
      }
    }
  }
}

void fm_refiner::update_gains_0_to_1(int v) {
#ifdef DEBUG_FM_REFINER
  assert(v >= 0 && v < numVertices);
  assert(partitionVector[v] == 0);
#endif

  int i;
  int j;
  int adjVertex;

  int hEdge;
  int endOffset;
  int destOffset;
  int sourceOffset;
  int endHedgeOffset;

  part_weights_[0] -= vWeight[v];
  part_weights_[1] += vWeight[v];

  partition_vector_[v] = 1;

  endOffset = vOffsets[v + 1];

  for (i = vOffsets[v]; i < endOffset; ++i) {
    hEdge = vToHedges[i];
    sourceOffset = Shiftl(hEdge, 1);
    destOffset = Or(Shiftl(hEdge, 1), 0x1);

    // ###
    // now if the hyperedge has no vertices in
    // destination part, increment the gains of
    // other free vertices in the hyperedge
    // ###

    if (vertex_in_part_[destOffset] == 0) {
      endHedgeOffset = hEdgeOffsets[hEdge + 1];

      for (j = hEdgeOffsets[hEdge]; j < endHedgeOffset; ++j) {
        adjVertex = pinList[j];

#ifdef DEBUG_FM_REFINER
        assert(adjVertex >= 0 && adjVertex < numVertices);
#endif
        if (adjVertex != v && partition_vector_[adjVertex] == 0) {
          vertex_gains_[adjVertex] += hEdgeWeight[hEdge];
          move_to_bucket_array(0, vertex_gains_[adjVertex], adjVertex);
        }
      }
    }

    // ###
    // move v from sP to dP for hyperedge
    // ###

    --vertex_in_part_[sourceOffset];
    ++vertex_in_part_[destOffset];

    // ###
    // if the hyperedge has one vertex in the sP
    // increment the gain of this vertex if free
    // ###

    if (vertex_in_part_[sourceOffset] == 1) {
      endHedgeOffset = hEdgeOffsets[hEdge + 1];

      for (j = hEdgeOffsets[hEdge]; j < endHedgeOffset; ++j) {
        adjVertex = pinList[j];

#ifdef DEBUG_REFINER
        assert(adjVertex >= 0 && adjVertex < numVertices);
#endif
        if (adjVertex != v && partition_vector_[adjVertex] == 0) {
          vertex_gains_[adjVertex] += hEdgeWeight[hEdge];
          move_to_bucket_array(0, vertex_gains_[adjVertex], adjVertex);
          break;
        }
      }
    }
  }
}

void fm_refiner::undo_move(int v) {
  int vPart = partition_vector_[v];
  int hEdge;

  int i;
  int endOffset = vOffsets[v + 1];
  int hEdgeOffset;
  int othPart = Xor(vPart, 0x1);

#ifdef DEBUG_FM_REFINER
  assert(v < numVertices && v >= 0);
  assert(vPart == 0 || vPart == 1);
#endif

  part_weights_[vPart] -= vWeight[v];
  part_weights_[othPart] += vWeight[v];

  partition_vector_[v] = othPart;

  for (i = vOffsets[v]; i < endOffset; ++i) {
    hEdge = vToHedges[i];
    hEdgeOffset = Shiftl(hEdge, 1);

    --vertex_in_part_[Or(hEdgeOffset, vPart)];
    ++vertex_in_part_[Or(hEdgeOffset, othPart)];
  }
}

void fm_refiner::refine(serial::hypergraph &h) {
  int totalGain;
  int gain;
  int i;
  int bestCutsize = LARGE_CONSTANT;
  int maxCutsize = 0;

  load_for_refinement(h);
  set_ee_threshold();
  build_buckets();

  for (i = 0; i < numPartitions; ++i) {
    set_partition_vector(i);
    initialise_partition_structure();

    totalGain = 0;

    if (part_weights_[0] > maximum_part_weight_)
      totalGain += rebalancing_pass(0);

    if (part_weights_[1] > maximum_part_weight_)
      totalGain += rebalancing_pass(1);

    if (part_weights_[0] <= maximum_part_weight_ && part_weights_[1] <= maximum_part_weight_) {
      do {
        gain = fm_pass();
        totalGain += gain;
      } while (gain > 0);
    }

    partitionCutsizes[i] -= totalGain;

    if (partitionCutsizes[i] < bestCutsize)
      bestCutsize = partitionCutsizes[i];

    if (partitionCutsizes[i] > maxCutsize)
      maxCutsize = partitionCutsizes[i];

    restore_buckets();
  }

  destroy_buckets();
}

int fm_refiner::fm_pass() {
  prepare_vertex_gains();
  loaded_.unset();

  int vGain;
  int numMoves = 0;
  int bestVertex;

  int i;
  int numOfBestMove = 0;
  int movesSincePosGain = 0;
  int gainSum = 0;
  int gainBest = 0;
  int bestVertexDP;

  while (numMoves < numVertices) {
    bestVertex = choose_maximum_gain_vertex();

    if (bestVertex != -1) {
      vGain = vertex_gains_[bestVertex];
      gainSum += vGain;
      bestVertexDP = partition_vector_[bestVertex];

      if (vGain <= 0)
        ++movesSincePosGain;

#ifdef DEBUG_FM_REFINER
      for (int ij = 0; ij < 2; ++ij)
        if (numBucketsInArray[ij] > 0) {
          int index = maxGainEntries[ij];
          assert((*bucketArrays[ij])[index].next);
        }
#endif
      remove_from_bucket_array(bestVertex, bestVertexDP, vGain);
      update_gains(bestVertex);
      loaded_.set(bestVertex);
      move_list_[numMoves++] = bestVertex;

      if (gainSum > gainBest) {
        gainBest = gainSum;
        numOfBestMove = numMoves;
      }

      if (movesSincePosGain > maximum_non_positive_moves_)
        break;
    } else {
      // ###
      // # no more 'feasible' vertex moves
      // ###

      break;
    }
  }

  // ###
  // # take back moves from best partial sum onwards
  // ###

  for (i = numMoves - 1; i >= numOfBestMove; --i)
    undo_move(move_list_[i]);

  remove_unlocked_from_bucket_arrays();

  return gainBest;
}

int fm_refiner::rebalancing_pass(int largePart) {
  int gain = 0;
  int vertexGain;
  int bestVertex;

  if (And(largePart, 0x1)) {
    initialise_gains_1_to_0();

    while (part_weights_[1] > maximum_part_weight_) {
      bestVertex = choose_legal_move(0x1);

      if (bestVertex == -1)
        break;

      vertexGain = vertex_gains_[bestVertex];
      gain += vertexGain;

      remove_from_bucket_array(bestVertex, 0x1, vertexGain);
      update_gains_1_to_0(bestVertex);
    }

    remove_buckets_from_1();
  } else {
    initialise_gains_0_to_1();

    while (part_weights_[0] > maximum_part_weight_) {
      bestVertex = choose_legal_move(0);

      if (bestVertex == -1)
        break;

      vertexGain = vertex_gains_[bestVertex];
      gain += vertexGain;

      remove_from_bucket_array(bestVertex, 0, vertexGain);
      update_gains_0_to_1(bestVertex);
    }

    remove_buckets_from_0();
  }

  return gain;
}

int fm_refiner::choose_maximum_gain_vertex() {
  int _0to1V = -1;
  int _1to0V = -1;
  int _0to1Gain = -maximum_possible_gain_;
  int _1to0Gain = -maximum_possible_gain_;

  int index;
  int v;

  bucket *bucket;

  if (number_of_buckets_in_array_[0] > 0) {
    index = maximum_gain_entries_[0];

#ifdef DEBUG_FM_REFINER
    assert((*bucketArrays[0])[index].next);
#endif

    bucket = (*bucket_arrays_[0])[index].next;
    v = bucket->vertex_id;

    while (part_weights_[1] + vWeight[v] > maximum_part_weight_) {
      bucket = bucket->next;

      if (!bucket)
        break;
      else
        v = bucket->vertex_id;
    }

    if (bucket) {
#ifdef DEBUG_FM_REFINER
      assert(v != -1);
#endif
      _0to1Gain = maximum_gains_[0];
      _0to1V = v;
    }
  }

  if (number_of_buckets_in_array_[1] > 0) {
    index = maximum_gain_entries_[1];

#ifdef DEBUG_FM_REFINER
    assert((*bucketArrays[1])[index].next);
#endif

    bucket = (*bucket_arrays_[1])[index].next;
    v = bucket->vertex_id;

    while (part_weights_[0] + vWeight[v] > maximum_part_weight_) {
      bucket = bucket->next;

      if (!bucket)
        break;
      else
        v = bucket->vertex_id;
    }

    if (bucket) {
#ifdef DEBUG_FM_REFINER
      assert(v != -1);
#endif
      _1to0Gain = maximum_gains_[1];
      _1to0V = v;
    }
  }

  if (_0to1V != -1 && _1to0V != -1) {
    if (_0to1Gain > _1to0Gain)
      return _0to1V;
    else if (_1to0Gain > _0to1Gain)
      return _1to0V;
    else if (part_weights_[0] + vWeight[_1to0V] >
             part_weights_[1] + vWeight[_0to1V])
      return _0to1V;
    else
      return _1to0V;
  } else {
    if (_0to1V != -1) {
      return _0to1V;
    } else {
      if (_1to0V != -1) {
        return _1to0V;
      } else {
        return -1;
      }
    }
  }
}

int fm_refiner::choose_legal_move(int sP) {
#ifdef DEBUG_FM_REFINER
  assert(numBucketsInArray[sP] > 0);
#endif

  int v;
  int index;

  bucket *b;

  if (And(sP, 0x1)) {
    index = maximum_gain_entries_[1];

#ifdef DEBUG_FM_REFINER
    assert((*bucketArrays[1])[index].next);
#endif

    b = (*bucket_arrays_[1])[index].next;
    v = b->vertex_id;

#ifdef DEBUG_FM_REFINER
    assert(v >= 0 && v < numVertices);
#endif

    while (part_weights_[0] + vWeight[v] > maximum_part_weight_) {
      b = b->next;

      while (!b && index > 0) {
        --index;
#ifdef DEBUG_FM_REFINER
        assert(index >= 0);
#endif
        b = (*bucket_arrays_[1])[index].next;
      }

      if (!b)
        return -1;

      v = b->vertex_id;

#ifdef DEBUG_FM_REFINER
      assert(v >= 0 && v < numVertices);
#endif
    }
  } else {
    index = maximum_gain_entries_[0];

#ifdef DEBUG_FM_REFINER
    assert((*bucketArrays[0])[index].next);
#endif

    b = (*bucket_arrays_[0])[index].next;
    v = b->vertex_id;

#ifdef DEBUG_FM_REFINER
    assert(v >= 0 && v < numVertices);
#endif

    while (part_weights_[1] + vWeight[v] > maximum_part_weight_) {
      b = b->next;

      while (!b && index > 0) {
        --index;
#ifdef DEBUG_FM_REFINER
        assert(index >= 0);
#endif
        b = (*bucket_arrays_[0])[index].next;
      }

      if (!b)
        return -1;

      v = b->vertex_id;

#ifdef DEBUG_FM_REFINER
      assert(v >= 0 && v < numVertices);
#endif
    }
  }

  return v;
}

}  // namespace serial
}  // namespace parkway
