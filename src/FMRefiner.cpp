#ifndef _FM_REFINER_CPP
#define _FM_REFINER_CPP

// ### FMRefiner.cpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 30/11/2004: Last Modified
//
// ###

#include "FMRefiner.hpp"

FMRefiner::FMRefiner(int max, int insMethod, int ee, int dL) : Refiner(dL) {
  bucketArraysLen = 0;
  maxPossGain = 0;
  maxNonPosMoves = 0;
  buckets.reserve(0);
  moveList.reserve(0);
  vertexGains.reserve(0);
  vInPart.reserve(0);
  locked.reserve(0);

  eeThreshold = ee;
  insertMethod = insMethod;
  maxPartWt = max;
  numParts = 2;
  partWeights.reserve(2);
  bucketArrays.reserve(2);
  numBucketsInArray.reserve(2);
  maxGainEntries.reserve(2);
  maxGains.reserve(2);

  bucketArrays[0] = nullptr;
  bucketArrays[1] = nullptr;
}

FMRefiner::~FMRefiner() {
  DynaMem::deletePtr<NodeArray>(bucketArrays[0]);
  DynaMem::deletePtr<NodeArray>(bucketArrays[1]);
}

void FMRefiner::dispRefinerOptions(std::ostream &out) const {
  switch (dispOption) {
  case SILENT:
    break;

  default:

    out << "|- FM:"
        << " qDis = ";
    printQdis(out);
    out << " eeT = " << eeThreshold << " acc = " << acceptProp << std::endl
        << "|" << std::endl;
    break;
  }
}

void FMRefiner::printQdis(std::ostream &out) const {
  switch (insertMethod) {
  case FIFO:
    out << "FIFO";
    break;

  default:
    out << "LIFO";
    break;
  }
}

void FMRefiner::buildBuckets() {
  int i;
  int j;

  int maxVertDeg = 0;
  int maxHedgeWt = 0;

  buckets.reserve(numVertices);
  vertexGains.reserve(numVertices);
  locked.reserve(numVertices);
  moveList.reserve(numVertices);
  vInPart.reserve(Shiftl(numHedges, 1));

  for (i = 0; i < numVertices; ++i) {
    buckets[i] = new bucket_node;
    buckets[i]->vertex_id = i;
    buckets[i]->previous = nullptr;
    buckets[i]->next = nullptr;

    j = vOffsets[i + 1] - vOffsets[i];

    if (j > maxVertDeg)
      maxVertDeg = j;
  }

  for (i = 0; i < numHedges; ++i) {
    if (hEdgeWeight[i] > maxHedgeWt) {
      maxHedgeWt = hEdgeWeight[i];
    }
  }

  maxPossGain = maxVertDeg * maxHedgeWt;
  bucketArraysLen = Or(Shiftl(maxPossGain, 1), 0x1);

  bucketArrays[0] = new NodeArray(bucketArraysLen);
  bucketArrays[1] = new NodeArray(bucketArraysLen);

  for (i = 0; i < bucketArraysLen; ++i) {
    (*bucketArrays[0])[i].next = nullptr;
    (*bucketArrays[0])[i].vertex_id = -1;
    (*bucketArrays[1])[i].next = nullptr;
    (*bucketArrays[1])[i].vertex_id = -1;
  }
}

void FMRefiner::restoreBuckets() {
  int i;

  for (i = 0; i < bucketArraysLen; ++i) {
    (*bucketArrays[0])[i].next = nullptr;
    (*bucketArrays[0])[i].vertex_id = -1;
    (*bucketArrays[1])[i].next = nullptr;
    (*bucketArrays[1])[i].vertex_id = -1;
  }

  for (i = 0; i < numVertices; ++i) {
    buckets[i]->previous = nullptr;
    buckets[i]->next = nullptr;
  }
}

void FMRefiner::destroyBuckets() {
  int i;

  for (i = 0; i < numVertices; ++i)
    DynaMem::deletePtr<bucket_node>(buckets[i]);

  for (i = 0; i < 2; ++i)
    DynaMem::deletePtr<NodeArray>(bucketArrays[i]);
}

void FMRefiner::initPartitionStruct() {
  int i;
  int j;
  int hEdgeOffset;

  int endOffset;

  partWeights[0] = 0;
  partWeights[1] = 0;

  for (i = 0; i < numVertices; ++i) {
    partWeights[partitionVector[i]] += vWeight[i];
  }

  for (i = 0; i < numHedges; ++i) {
    hEdgeOffset = Shiftl(i, 1);

    vInPart[hEdgeOffset] = 0;
    vInPart[Or(hEdgeOffset, 0x1)] = 0;

    endOffset = hEdgeOffsets[i + 1];

    for (j = hEdgeOffsets[i]; j < endOffset; ++j) {
      ++vInPart[Or(hEdgeOffset, partitionVector[pinList[j]])];
    }
  }
}

void FMRefiner::prepVertexGains() {
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
    maxGains[i] = -maxPossGain;
    maxGainEntries[i] = 0;
    numBucketsInArray[i] = 0;
  }

  for (i = 0; i < numVertices; ++i) {
    endOffset = vOffsets[i + 1];
    vPart = partitionVector[i];

#ifdef DEBUG_FM_REFINER
    assert(vPart == 0 || vPart == 1);
#endif

    vertexGains[i] = 0;

    for (j = vOffsets[i]; j < endOffset; ++j) {
      hEdge = vToHedges[j];
      hEdgeOffset = Shiftl(hEdge, 1);

      if (vInPart[Or(hEdgeOffset, vPart)] == 1)
        vertexGains[i] += hEdgeWeight[hEdge];

      if (vInPart[Or(hEdgeOffset, Xor(vPart, 1))] == 0)
        vertexGains[i] -= hEdgeWeight[hEdge];
    }

    moveToBucketArray(vPart, vertexGains[i], i);
    ++numBucketsInArray[vPart];
  }
}

void FMRefiner::initGains1to0() {
  int i;
  int j;

  int hEdge;
  int endOffset;

  // ###
  // compute vertex gains in 1 to 0 direction
  // ###

  for (i = 0; i < 2; ++i) {
    maxGains[i] = -maxPossGain;
    maxGainEntries[i] = 0;
    numBucketsInArray[i] = 0;
  }

  for (i = 0; i < numVertices; ++i) {
    if (partitionVector[i] == 1) {
      endOffset = vOffsets[i + 1];
      vertexGains[i] = 0;

      for (j = vOffsets[i]; j < endOffset; ++j) {
        hEdge = vToHedges[j];

        if (vInPart[Shiftl(hEdge, 1)] == 0)
          vertexGains[i] -= hEdgeWeight[hEdge];

        if (vInPart[Or(Shiftl(hEdge, 1), 0x1)] == 1)
          vertexGains[i] += hEdgeWeight[hEdge];
      }

      moveToBucketArray(1, vertexGains[i], i);

#ifdef DEBUG_FM_REFINER
      assert(buckets[i]->prev != nullptr);
#endif
      ++numBucketsInArray[1];
    }
  }
}

void FMRefiner::initGains0to1() {
  int i;
  int j;

  int hEdge;
  int endOffset;

  // ###
  // compute vertex gains in 0 to 1 direction
  // ###

  for (i = 0; i < 2; ++i) {
    maxGains[i] = -maxPossGain;
    maxGainEntries[i] = 0;
    numBucketsInArray[i] = 0;
  }

  for (i = 0; i < numVertices; ++i) {
    if (partitionVector[i] == 0) {
      endOffset = vOffsets[i + 1];
      vertexGains[i] = 0;

      for (j = vOffsets[i]; j < endOffset; ++j) {
        hEdge = vToHedges[j];

        if (vInPart[Or(Shiftl(hEdge, 1), 0x1)] == 0)
          vertexGains[i] -= hEdgeWeight[hEdge];

        if (vInPart[Shiftl(hEdge, 1)] == 1)
          vertexGains[i] += hEdgeWeight[hEdge];
      }

      moveToBucketArray(0, vertexGains[i], i);

#ifdef DEBUG_FM_REFINER
      assert(buckets[i]->prev != nullptr);
#endif
      ++numBucketsInArray[0];
    }
  }
}

void FMRefiner::removeBucketsFrom1() {
  int i;
  bucket_node *b;

  bucket_node *array = bucketArrays[1]->data();

  int bucketArrayLen = Or(Shiftl(maxPossGain, 1), 0x1);

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

  numBucketsInArray[1] = 0;
}

void FMRefiner::removeBucketsFrom0() {
  int i;
  bucket_node *b;

  bucket_node *array = bucketArrays[0]->data();

  int bucketArrayLen = Or(Shiftl(maxPossGain, 1), 0x1);

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

  numBucketsInArray[0] = 0;
}

void FMRefiner::removeUnlockedFromBucketArrays() {
  int i;
  bucket_node *b;

  for (i = 0; i < numVertices; ++i) {
    if (!locked(i)) {
      b = buckets[i];

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

  numBucketsInArray[0] = 0;
  numBucketsInArray[1] = 0;
}

void FMRefiner::moveToBucketArray(int vPart, int vGain, int v) {
#ifdef DEBUGH_REFINER
  assert(vPart == 0 || vPart == 1);
  assert(v >= 0 && v < numVertices);
#endif

  int arrayIndex = vGain + maxPossGain;

  bucket_node *b = buckets[v];
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

  if (vGain > maxGains[vPart]) {
    maxGainEntries[vPart] = arrayIndex;
    maxGains[vPart] = vGain;
  }

  bucket = &(*bucketArrays[vPart])[arrayIndex];

  switch (insertMethod) {
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

  if (maxGains[vPart] > vGain) {
    int index = maxGainEntries[vPart];

    while ((*bucketArrays[vPart])[index].next == nullptr) {
#ifdef DEBUG_FM_REFINER
      assert(index >= 0);
#endif
      --index;
    }

    maxGainEntries[vPart] = index;
    maxGains[vPart] = index - maxPossGain;

#ifdef DEBUG_FM_REFINER
    assert(index >= 0);
    assert((*bucketArrays[vPart])[index].next);
#endif
  }
}

void FMRefiner::removeFromBucketArray(int v, int vPart, int vGain) {
#ifdef DEBUG_FM_REFINER
  assert(vPart == 0 || vPart == 1);
  assert(v >= 0 && v < numVertices);
#endif

  bucket_node *b = buckets[v];

#ifdef DEBUG_FM_REFINER
  assert(b->prev);
#endif

  b->previous->next = b->next;

  if (b->next)
    b->next->previous = b->previous;

  b->previous = nullptr;
  b->next = nullptr;

  --numBucketsInArray[vPart];

#ifdef DEBUG_FM_REFINER
  assert(numBucketsInArray[vPart] >= 0);
#endif

  if (maxGains[vPart] == vGain) {
    int index = vGain + maxPossGain;

    while ((*bucketArrays[vPart])[index].next == nullptr && index > 0)
      --index;

    maxGainEntries[vPart] = index;
    maxGains[vPart] = index - maxPossGain;

#ifdef DEBUG_FM_REFINER
    assert(index >= 0);
#endif
  }
}

void FMRefiner::adjustMaxGainPtr(int dP) {
  if (numBucketsInArray[dP] > 0) {
    int index = Shiftl(maxPossGain, 1);

    while ((*bucketArrays[dP])[index].next == nullptr) {
#ifdef DEBUG_FM_REFINER
      assert(index > 0);
#endif
      index--;
    }

    maxGainEntries[dP] = index;
    maxGains[dP] = index - maxPossGain;

#ifdef DEBUG_FM_REFINER
    assert(index >= 0);
    assert((*bucketArrays[dP])[index].next);
#endif
  }
}

void FMRefiner::updateGains(int v) {
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
  int sP = partitionVector[v];
  int dP = sP ^ 0x1;

  partWeights[sP] -= vWeight[v];
  partWeights[dP] += vWeight[v];

  partitionVector[v] = dP;

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

    if (vInPart[destOffset] == 0) {
      endHedgeOffset = hEdgeOffsets[hEdge + 1];

      for (j = hEdgeOffsets[hEdge]; j < endHedgeOffset; ++j) {
        adjVertex = pinList[j];

#ifdef DEBUG_FM_REFINER
        if (adjVertex != v)
          assert(partitionVector[adjVertex] == sP);
#endif
        if (adjVertex != v && !locked(adjVertex)) {
          vertexGains[adjVertex] += hEdgeWeight[hEdge];
          moveToBucketArray(sP, vertexGains[adjVertex], adjVertex);
        }
      }
    }

    // ###
    // if the hyperedge has one vertex in the
    // destination part and it is free, decrement
    // its gain
    // ###

    if (vInPart[destOffset] == 1) {
      endHedgeOffset = hEdgeOffsets[hEdge + 1];

      for (j = hEdgeOffsets[hEdge]; j < endHedgeOffset; ++j) {
        adjVertex = pinList[j];

        if (adjVertex != v && partitionVector[adjVertex] == dP &&
            !locked(adjVertex)) {
          vertexGains[adjVertex] -= hEdgeWeight[hEdge];
          moveToBucketArray(dP, vertexGains[adjVertex], adjVertex);
          break;
        }
      }
    }

    // ###
    // move v from sP to dP for hyperedge
    // ###

    --vInPart[sourceOffset];
    ++vInPart[destOffset];

    // ###
    // if the hyperedge has no vertices in the sP
    // decrement gains of all free vertices
    // ###

    if (vInPart[sourceOffset] == 0) {
      endHedgeOffset = hEdgeOffsets[hEdge + 1];

      for (j = hEdgeOffsets[hEdge]; j < endHedgeOffset; ++j) {
        adjVertex = pinList[j];

#ifdef DEBUG_FM_REFINER
        if (adjVertex != v)
          assert(partitionVector[adjVertex] == dP);
#endif

        if (adjVertex != v && !locked(adjVertex)) {
          vertexGains[adjVertex] -= hEdgeWeight[hEdge];
          moveToBucketArray(dP, vertexGains[adjVertex], adjVertex);
        }
      }
    }

    // ###
    // if the hyperedge has one vertex in the sP
    // increment the gain of this vertex if free
    // ###

    if (vInPart[sourceOffset] == 1) {
      endHedgeOffset = hEdgeOffsets[hEdge + 1];

      for (j = hEdgeOffsets[hEdge]; j < endHedgeOffset; ++j) {
        adjVertex = pinList[j];

        if (adjVertex != v && partitionVector[adjVertex] == sP &&
            !locked(adjVertex)) {
          vertexGains[adjVertex] += hEdgeWeight[hEdge];
          moveToBucketArray(sP, vertexGains[adjVertex], adjVertex);
          break;
        }
      }
    }
  }
}

void FMRefiner::updateGains1_0(int v) {
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

  partWeights[1] -= vWeight[v];
  partWeights[0] += vWeight[v];

  partitionVector[v] = 0;

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

    if (vInPart[destOffset] == 0) {
      endHedgeOffset = hEdgeOffsets[hEdge + 1];

      for (j = hEdgeOffsets[hEdge]; j < endHedgeOffset; ++j) {
        adjVertex = pinList[j];

#ifdef DEBUG_FM_REFINER
        assert(adjVertex >= 0 && adjVertex < numVertices);
#endif
        if (adjVertex != v && partitionVector[adjVertex] == 1) {
          vertexGains[adjVertex] += hEdgeWeight[hEdge];
          moveToBucketArray(1, vertexGains[adjVertex], adjVertex);
        }
      }
    }

    // ###
    // move v from sP to dP for hyperedge
    // ###

    --vInPart[sourceOffset];
    ++vInPart[destOffset];

    // ###
    // if the hyperedge has one vertex in the sP
    // increment the gain of this vertex if free
    // ###

    if (vInPart[sourceOffset] == 1) {
      endHedgeOffset = hEdgeOffsets[hEdge + 1];

      for (j = hEdgeOffsets[hEdge]; j < endHedgeOffset; ++j) {
        adjVertex = pinList[j];

#ifdef DEBUG_REFINER
        assert(adjVertex >= 0 && adjVertex < numVertices);
#endif
        if (adjVertex != v && partitionVector[adjVertex] == 1) {
          vertexGains[adjVertex] += hEdgeWeight[hEdge];
          moveToBucketArray(1, vertexGains[adjVertex], adjVertex);
          break;
        }
      }
    }
  }
}

void FMRefiner::updateGains0_1(int v) {
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

  partWeights[0] -= vWeight[v];
  partWeights[1] += vWeight[v];

  partitionVector[v] = 1;

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

    if (vInPart[destOffset] == 0) {
      endHedgeOffset = hEdgeOffsets[hEdge + 1];

      for (j = hEdgeOffsets[hEdge]; j < endHedgeOffset; ++j) {
        adjVertex = pinList[j];

#ifdef DEBUG_FM_REFINER
        assert(adjVertex >= 0 && adjVertex < numVertices);
#endif
        if (adjVertex != v && partitionVector[adjVertex] == 0) {
          vertexGains[adjVertex] += hEdgeWeight[hEdge];
          moveToBucketArray(0, vertexGains[adjVertex], adjVertex);
        }
      }
    }

    // ###
    // move v from sP to dP for hyperedge
    // ###

    --vInPart[sourceOffset];
    ++vInPart[destOffset];

    // ###
    // if the hyperedge has one vertex in the sP
    // increment the gain of this vertex if free
    // ###

    if (vInPart[sourceOffset] == 1) {
      endHedgeOffset = hEdgeOffsets[hEdge + 1];

      for (j = hEdgeOffsets[hEdge]; j < endHedgeOffset; ++j) {
        adjVertex = pinList[j];

#ifdef DEBUG_REFINER
        assert(adjVertex >= 0 && adjVertex < numVertices);
#endif
        if (adjVertex != v && partitionVector[adjVertex] == 0) {
          vertexGains[adjVertex] += hEdgeWeight[hEdge];
          moveToBucketArray(0, vertexGains[adjVertex], adjVertex);
          break;
        }
      }
    }
  }
}

void FMRefiner::unmakeMove(int v) {
  int vPart = partitionVector[v];
  int hEdge;

  int i;
  int endOffset = vOffsets[v + 1];
  int hEdgeOffset;
  int othPart = Xor(vPart, 0x1);

#ifdef DEBUG_FM_REFINER
  assert(v < numVertices && v >= 0);
  assert(vPart == 0 || vPart == 1);
#endif

  partWeights[vPart] -= vWeight[v];
  partWeights[othPart] += vWeight[v];

  partitionVector[v] = othPart;

  for (i = vOffsets[v]; i < endOffset; ++i) {
    hEdge = vToHedges[i];
    hEdgeOffset = Shiftl(hEdge, 1);

    --vInPart[Or(hEdgeOffset, vPart)];
    ++vInPart[Or(hEdgeOffset, othPart)];
  }
}

void FMRefiner::refine(serial::hypergraph &h) {
  int totalGain;
  int gain;
  int i;
  int bestCutsize = LARGE_CONSTANT;
  int maxCutsize = 0;

  load_for_refinement(h);
  setEEThreshold();
  buildBuckets();

  for (i = 0; i < numPartitions; ++i) {
    setPartitionVector(i);
    initPartitionStruct();

    totalGain = 0;

    if (partWeights[0] > maxPartWt)
      totalGain += doRebalancingPass(0);

    if (partWeights[1] > maxPartWt)
      totalGain += doRebalancingPass(1);

    if (partWeights[0] <= maxPartWt && partWeights[1] <= maxPartWt) {
      do {
        gain = doFidMatPass();
        totalGain += gain;
      } while (gain > 0);
    }

    partitionCutsizes[i] -= totalGain;

    if (partitionCutsizes[i] < bestCutsize)
      bestCutsize = partitionCutsizes[i];

    if (partitionCutsizes[i] > maxCutsize)
      maxCutsize = partitionCutsizes[i];

    restoreBuckets();
  }

  destroyBuckets();
}

int FMRefiner::doFidMatPass() {
  prepVertexGains();
  locked.unset();

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
    bestVertex = chooseMaxGainVertex();

    if (bestVertex != -1) {
      vGain = vertexGains[bestVertex];
      gainSum += vGain;
      bestVertexDP = partitionVector[bestVertex];

      if (vGain <= 0)
        ++movesSincePosGain;

#ifdef DEBUG_FM_REFINER
      for (int ij = 0; ij < 2; ++ij)
        if (numBucketsInArray[ij] > 0) {
          int index = maxGainEntries[ij];
          assert((*bucketArrays[ij])[index].next);
        }
#endif
      removeFromBucketArray(bestVertex, bestVertexDP, vGain);
      updateGains(bestVertex);
      locked.set(bestVertex);
      moveList[numMoves++] = bestVertex;

      if (gainSum > gainBest) {
        gainBest = gainSum;
        numOfBestMove = numMoves;
      }

      if (movesSincePosGain > maxNonPosMoves)
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
    unmakeMove(moveList[i]);

  removeUnlockedFromBucketArrays();

  return gainBest;
}

int FMRefiner::doRebalancingPass(int largePart) {
  int gain = 0;
  int vertexGain;
  int bestVertex;

  if (And(largePart, 0x1)) {
    initGains1to0();

    while (partWeights[1] > maxPartWt) {
      bestVertex = chooseLegalMove(0x1);

      if (bestVertex == -1)
        break;

      vertexGain = vertexGains[bestVertex];
      gain += vertexGain;

      removeFromBucketArray(bestVertex, 0x1, vertexGain);
      updateGains1_0(bestVertex);
    }

    removeBucketsFrom1();
  } else {
    initGains0to1();

    while (partWeights[0] > maxPartWt) {
      bestVertex = chooseLegalMove(0);

      if (bestVertex == -1)
        break;

      vertexGain = vertexGains[bestVertex];
      gain += vertexGain;

      removeFromBucketArray(bestVertex, 0, vertexGain);
      updateGains0_1(bestVertex);
    }

    removeBucketsFrom0();
  }

  return gain;
}

int FMRefiner::chooseMaxGainVertex() {
  int _0to1V = -1;
  int _1to0V = -1;
  int _0to1Gain = -maxPossGain;
  int _1to0Gain = -maxPossGain;

  int index;
  int v;

  bucket *bucket;

  if (numBucketsInArray[0] > 0) {
    index = maxGainEntries[0];

#ifdef DEBUG_FM_REFINER
    assert((*bucketArrays[0])[index].next);
#endif

    bucket = (*bucketArrays[0])[index].next;
    v = bucket->vertex_id;

    while (partWeights[1] + vWeight[v] > maxPartWt) {
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
      _0to1Gain = maxGains[0];
      _0to1V = v;
    }
  }

  if (numBucketsInArray[1] > 0) {
    index = maxGainEntries[1];

#ifdef DEBUG_FM_REFINER
    assert((*bucketArrays[1])[index].next);
#endif

    bucket = (*bucketArrays[1])[index].next;
    v = bucket->vertex_id;

    while (partWeights[0] + vWeight[v] > maxPartWt) {
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
      _1to0Gain = maxGains[1];
      _1to0V = v;
    }
  }

  if (_0to1V != -1 && _1to0V != -1) {
    if (_0to1Gain > _1to0Gain)
      return _0to1V;
    else if (_1to0Gain > _0to1Gain)
      return _1to0V;
    else if (partWeights[0] + vWeight[_1to0V] >
             partWeights[1] + vWeight[_0to1V])
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

int FMRefiner::chooseLegalMove(int sP) {
#ifdef DEBUG_FM_REFINER
  assert(numBucketsInArray[sP] > 0);
#endif

  int v;
  int index;

  bucket *b;

  if (And(sP, 0x1)) {
    index = maxGainEntries[1];

#ifdef DEBUG_FM_REFINER
    assert((*bucketArrays[1])[index].next);
#endif

    b = (*bucketArrays[1])[index].next;
    v = b->vertex_id;

#ifdef DEBUG_FM_REFINER
    assert(v >= 0 && v < numVertices);
#endif

    while (partWeights[0] + vWeight[v] > maxPartWt) {
      b = b->next;

      while (!b && index > 0) {
        --index;
#ifdef DEBUG_FM_REFINER
        assert(index >= 0);
#endif
        b = (*bucketArrays[1])[index].next;
      }

      if (!b)
        return -1;

      v = b->vertex_id;

#ifdef DEBUG_FM_REFINER
      assert(v >= 0 && v < numVertices);
#endif
    }
  } else {
    index = maxGainEntries[0];

#ifdef DEBUG_FM_REFINER
    assert((*bucketArrays[0])[index].next);
#endif

    b = (*bucketArrays[0])[index].next;
    v = b->vertex_id;

#ifdef DEBUG_FM_REFINER
    assert(v >= 0 && v < numVertices);
#endif

    while (partWeights[1] + vWeight[v] > maxPartWt) {
      b = b->next;

      while (!b && index > 0) {
        --index;
#ifdef DEBUG_FM_REFINER
        assert(index >= 0);
#endif
        b = (*bucketArrays[0])[index].next;
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

#endif
