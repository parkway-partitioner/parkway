
#ifndef _KWAY_GREEDY_REFINER_CPP
#define _KWAY_GREEDY_REFINER_CPP

// ### GreedyKwayRefiner.cpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 30/11/2004: Last Modified
//
// ###

#include "GreedyKwayRefiner.hpp"

GreedyKwayRefiner::GreedyKwayRefiner(int max, int nparts, double ave,
                                     double lim, int dL)
    : Refiner(dL) {
  maxPartWt = max;
  numParts = nparts;
  avePartWt = ave;
  limit = lim;
  numNonPosMoves = 0;
  partWeights.reserve(numParts);

  numNeighParts.reserve(0);
  neighboursOfV.reserve(0);
  neighboursOfVOffsets.reserve(0);
  hEdgeVinPart.reserve(0);
  hEdgeVinPartOffsets.reserve(0);
  vertices.reserve(0);
  vertSeen.reserve(0);
  seenVertices.reserve(0);
  partsSpanned.reserve(0);

#ifdef DEBUG_REFINER
  assert(limit >= 0 && limit <= 1.0);
#endif
}

GreedyKwayRefiner::~GreedyKwayRefiner() {}

void GreedyKwayRefiner::dispRefinerOptions(std::ostream &out) const {
  switch (dispOption) {
  case SILENT:
    break;

  default:

    out << "|- GKWAY:"
        << " lim = " << limit << std::endl
        << "|" << std::endl;
    break;
  }
}

void GreedyKwayRefiner::buildDataStructs() {
  int i;

  partsSpanned.reserve(numParts);
  vertices.reserve(numVertices);
  vertSeen.reserve(numVertices);
  seenVertices.reserve(numVertices);
  numNeighParts.reserve(numVertices);
  neighboursOfVOffsets.reserve(numVertices + 1);
  neighboursOfV.reserve(numVertices * numParts);

  hEdgeVinPartOffsets.reserve(numHedges + 1);
  hEdgeVinPart.reserve(numHedges * numParts);

  neighboursOfVOffsets[0] = 0;
  hEdgeVinPartOffsets[0] = 0;

  for (i = 1; i <= numVertices; ++i)
    neighboursOfVOffsets[i] = neighboursOfVOffsets[i - 1] + numParts;

  for (i = 1; i <= numHedges; ++i)
    hEdgeVinPartOffsets[i] = hEdgeVinPartOffsets[i - 1] + numParts;
}

void GreedyKwayRefiner::destroyDataStructs() {
  partsSpanned.reserve(0);
  vertices.reserve(0);
  vertSeen.reserve(0);
  seenVertices.reserve(0);
  numNeighParts.reserve(0);
  neighboursOfVOffsets.reserve(0);
  neighboursOfV.reserve(0);

  hEdgeVinPartOffsets.reserve(0);
  hEdgeVinPart.reserve(0);
}

int GreedyKwayRefiner::initDataStructs() {
  int i;
  int j;
  int ij;

  int v;
  int endIndex;
  int vertOffset;
  int partOffset;
  int numPartsSpanned;
  int vPart;
  int k_1cut = 0;

  for (i = 0; i < numParts; ++i)
    partWeights[i] = 0;

  for (i = 0; i < numVertices; ++i) {
#ifdef DEBUG_REFINER
    assert(partitionVector[i] >= 0 && partitionVector[i] < numParts);
#endif
    partWeights[partitionVector[i]] += vWeight[i];
    numNeighParts[i] = 0;
    vertSeen[i] = -1;
  }

  endIndex = neighboursOfVOffsets[numVertices];

  for (i = 0; i < endIndex; ++i)
    neighboursOfV[i] = 0;

  endIndex = hEdgeVinPartOffsets[numHedges];

  for (i = 0; i < endIndex; ++i)
    hEdgeVinPart[i] = 0;

  // ###
  // initialise the hyperedge structures
  // and the vertices structure
  // ###

  for (i = 0; i < numHedges; ++i) {
    endIndex = hEdgeOffsets[i + 1];
    numPartsSpanned = 0;

    // ###
    // update the hyperedge part distributions
    // ###

    for (j = hEdgeOffsets[i]; j < endIndex; ++j) {
      vPart = partitionVector[pinList[j]];
      partOffset = hEdgeVinPartOffsets[i] + vPart;

      if (hEdgeVinPart[partOffset] == 0) {
        partsSpanned[numPartsSpanned++] = vPart;
      }

      ++hEdgeVinPart[partOffset];
    }

    // ###
    // update the vertex neighbour structs
    // ###

    for (j = hEdgeOffsets[i]; j < endIndex; ++j) {
      v = pinList[j];
      vertOffset = neighboursOfVOffsets[v];

      for (ij = 0; ij < numPartsSpanned; ++ij) {
        vPart = partsSpanned[ij];

        if (neighboursOfV[vertOffset + vPart] == 0)
          ++numNeighParts[v];

        neighboursOfV[vertOffset + vPart] = 1;
      }
    }

    k_1cut += ((numPartsSpanned - 1) * hEdgeWeight[i]);
  }

  return k_1cut;
}

void GreedyKwayRefiner::updateAdjVertStat(int v, int sP, int bestDP) {
  int i;
  int j;
  int ij;
  int numVerticesSeen;
  int vertOffset;
  int vert;
  int hEdge;
  int hEdgeOffset;
  int neighOfVOffset;
  int othVOffset;
  int othHedge;
  int reduceNumNeighPartsOn;

  vertOffset = vOffsets[v + 1];
  numVerticesSeen = 0;

  partsSpanned.reserve(numParts);

  for (j = vOffsets[v]; j < vertOffset; ++j) {

    hEdge = vToHedges[j];
    hEdgeOffset = hEdgeOffsets[hEdge + 1];

    for (ij = hEdgeOffsets[hEdge]; ij < hEdgeOffset; ++ij) {
      vert = pinList[ij];

      if (vertSeen[vert] == -1 && vert != v) {
        neighOfVOffset = neighboursOfVOffsets[vert];

        if (neighboursOfV[neighOfVOffset + bestDP] == 0)
          ++numNeighParts[vert];

        neighboursOfV[neighOfVOffset + bestDP] = 1;

        if (partitionVector[vert] != sP) {
          if (neighboursOfV[neighOfVOffset + sP] > 0)
            reduceNumNeighPartsOn = 1;
          else
            reduceNumNeighPartsOn = 0;

          neighboursOfV[neighOfVOffset + sP] = 0;

          othVOffset = vOffsets[vert + 1];

          for (i = vOffsets[vert]; i < othVOffset; ++i) {
            othHedge = vToHedges[i];

            if (hEdgeVinPart[hEdgeVinPartOffsets[othHedge] + sP] > 0) {
              neighboursOfV[neighOfVOffset + sP] = 1;
              break;
            }
          }

          if (reduceNumNeighPartsOn && neighboursOfV[neighOfVOffset + sP] == 0)
            --numNeighParts[vert];
        }

        vertSeen[vert] = 1;
        seenVertices[numVerticesSeen] = vert;
        ++numVerticesSeen;
      }
    }
  }

  // ###
  // restore the 'seen' vertices structure
  // ###

  for (j = 0; j < numVerticesSeen; ++j)
    vertSeen[seenVertices[j]] = -1;
}

void GreedyKwayRefiner::refine(serial::hypergraph &h) {
  int totalGain = 0;
  int gain;
  int i;

  load_for_refinement(h);
  buildDataStructs();

  if (limit < 1.0) {
    setLimit();
  }

  for (i = 0; i < numPartitions; ++i) {
    setPartitionVector(i);
    initDataStructs();
    totalGain = 0;

    do {
      gain = runGreedyPass();
      totalGain += gain;
    } while (gain > 0);

    partitionCutsizes[i] -= totalGain;
  }

  destroyDataStructs();
}

void GreedyKwayRefiner::rebalance(serial::hypergraph &h) {
  int totalGain = 0;
  int gain;
  int i;

  load_for_refinement(h);
  buildDataStructs();

  for (i = 0; i < numPartitions; ++i) {
    setPartitionVector(i);
    initDataStructs();

    totalGain = runRebalancingPass();

    do {
      gain = runGreedyPass();
      totalGain += gain;
    } while (gain > 0);

    partitionCutsizes[i] -= totalGain;
  }

  destroyDataStructs();

#ifdef DEBUG_REFINER
  h.checkPartitions(numParts, maxPartWt);
#endif
}

int GreedyKwayRefiner::runGreedyPass() {
  int i;
  int j;
  int ij;
  int v;
  int sP;
  int tmp;
  int gain;
  int vGain;
  int posGain;
  int hEdge;
  int bestMove;
  int randomNum;
  int vertexWt;
  int vertOffset;
  int vNeighOffset;
  int hEdgeOffset;
  int neighOfVOffset;
  int numNonPos = 0;

  double currImbalance = 0;
  double bestImbalance;
  double posImbalance;

  for (i = 0; i < numVertices; ++i) {
    vertices[i] = i;
  }

  for (i = 0; i < numParts; ++i) {
    currImbalance += fabs(partWeights[i] - avePartWt);
  }

  i = numVertices;
  gain = 0;

  do {
    randomNum = RANDOM(0, i);
    v = vertices[randomNum];
    vertexWt = vWeight[v];
    sP = partitionVector[v];
    vGain = 0;
    bestImbalance = currImbalance;
    bestMove = -1;

    if (numNeighParts[v] > 1) {
      vNeighOffset = neighboursOfVOffsets[v];

      for (j = 0; j < numParts; ++j) {
        if (j != sP && neighboursOfV[vNeighOffset + j] > 0) {
          if (partWeights[j] + vertexWt <= maxPartWt) {
            posGain = 0;
            vertOffset = vOffsets[v + 1];

            for (ij = vOffsets[v]; ij < vertOffset; ++ij) {
              hEdge = vToHedges[ij];
              hEdgeOffset = hEdgeVinPartOffsets[hEdge];

              if (hEdgeVinPart[hEdgeOffset + sP] == 1)
                posGain += hEdgeWeight[hEdge];

              if (hEdgeVinPart[hEdgeOffset + j] == 0)
                posGain -= hEdgeWeight[hEdge];
            }

            posImbalance =
                currImbalance + fabs(partWeights[sP] - (vertexWt + avePartWt));
            posImbalance += fabs((partWeights[j] + vertexWt) - avePartWt);
            posImbalance -= fabs(partWeights[sP] - avePartWt);
            posImbalance -= fabs(partWeights[j] - avePartWt);

            if ((posGain > vGain) ||
                (posGain == vGain && posImbalance < bestImbalance)) {
              vGain = posGain;
              bestMove = j;
              bestImbalance = posImbalance;
            }
          }
        }
      }

      if (bestMove > -1) {
        vertOffset = vOffsets[v + 1];

        // ###
        // update the moved vertices' stats
        // ###

        neighOfVOffset = neighboursOfVOffsets[v];

        if (neighboursOfV[neighOfVOffset + bestMove] == 0)
          ++numNeighParts[v];

        neighboursOfV[neighOfVOffset + bestMove] = 1;
        neighboursOfV[neighOfVOffset + sP] = 0;

        for (j = vOffsets[v]; j < vertOffset; ++j) {
          // ###
          // update the hyperedge stats: (vInPart etc.)
          // ###

          hEdge = vToHedges[j];
          hEdgeOffset = hEdgeVinPartOffsets[hEdge];

          --hEdgeVinPart[hEdgeOffset + sP];
          ++hEdgeVinPart[hEdgeOffset + bestMove];

          if (hEdgeVinPart[hEdgeOffset + sP] > 0)
            neighboursOfV[neighOfVOffset + sP] = 1;
        }

        if (neighboursOfV[neighOfVOffset + sP] == 0)
          --numNeighParts[v];

        // ###
        // update the adj vertices stats:
        // (num neighbours in part etc.)
        // ###

        updateAdjVertStat(v, sP, bestMove);

        // ###
        // update other structs
        // ###

        partitionVector[v] = bestMove;
        partWeights[sP] -= vertexWt;
        partWeights[bestMove] += vertexWt;
        currImbalance = bestImbalance;

        // ###
        // finally, update the gain...
        // ###

        gain += vGain;

        if (vGain <= 0)
          ++numNonPos;
        else
          numNonPos = 0;

        if (limit < 1.0) {
          if (numNonPos > numNonPosMoves)
            break;
        }
      }
    }

    fswap(vertices[randomNum], vertices[i - 1], tmp);
    --i;
  } while (i > 0);

  return gain;
}

int GreedyKwayRefiner::runRebalancingPass() {
  int i;
  int j;

  int part;
  int numOverWeight = 0;

  dynamic_array<int> overWeight(numParts);

  for (i = 0; i < numParts; ++i) {
    if (partWeights[i] <= maxPartWt)
      overWeight[i] = 0;

    else {
      overWeight[i] = 1;
      ++numOverWeight;
    }
  }

  if (numOverWeight == 0)
    return 0;

  int vertex;
  int hEdge;
  int gain = 0;
  int vGain;
  int posGain;
  int bestMove;
  int vertexWt;
  int vertOffset;
  int hEdgeOffset;
  int neighOfVOffset;

  VNodePtr nodePtr;

  VNodePtrArray verticesInParts(numParts);
  VNodePtrArray vertexNodes(numVertices);

  for (i = 0; i < numVertices; ++i)
    vertices[i] = i;

  for (i = 0; i < numParts; ++i)
    verticesInParts[i] = nullptr;

  Funct::randomPermutation(vertices.data(), numVertices);

  for (i = 0; i < numVertices; ++i) {
    vertex = vertices[i];
    part = partitionVector[vertex];

    if (overWeight[part] == 0) {
      vertexNodes[vertex] = nullptr;
    } else {
      vertexNodes[vertex] = new VNode;
      vertexNodes[vertex]->vertexID = vertex;
      vertexNodes[vertex]->next = verticesInParts[part];
      verticesInParts[part] = vertexNodes[vertex];
    }
  }

  part = findHeaviestOverweight();

  while (part > -1) {
#ifdef DEBUG_REFINER
    assert(part >= 0 && part < numParts);
#endif

    nodePtr = verticesInParts[part];
    verticesInParts[part] = nodePtr->next;
    vertex = nodePtr->vertexID;
    vertexWt = vWeight[vertex];
    vGain = -LARGE_CONSTANT;
    bestMove = -1;

#ifdef DEBUG_REFINER
    assert(vertex >= 0 && vertex < numVertices);
    assert(partitionVector[vertex] == part);
#endif

    for (i = 0; i < numParts; ++i) {
      if (i != part && overWeight[i] == 0) {
        if (partWeights[i] + vertexWt <= maxPartWt) {
          posGain = 0;
          vertOffset = vOffsets[vertex + 1];

          for (j = vOffsets[vertex]; j < vertOffset; ++j) {
            hEdge = vToHedges[j];
            hEdgeOffset = hEdgeVinPartOffsets[hEdge];

            if (hEdgeVinPart[hEdgeOffset + part] == 1)
              posGain += hEdgeWeight[hEdge];

            if (hEdgeVinPart[hEdgeOffset + i] == 0)
              posGain -= hEdgeWeight[hEdge];
          }

          if (posGain > vGain) {
            vGain = posGain;
            bestMove = i;
          }
        }
      }
    }
#ifdef DEBUG_REFINER
    assert(bestMove != -1);
#endif

    vertOffset = vOffsets[vertex + 1];

    // ###
    // update the moved vertices' stats
    // ###

    neighOfVOffset = neighboursOfVOffsets[vertex];

    if (neighboursOfV[neighOfVOffset + bestMove] == 0)
      ++numNeighParts[vertex];

    neighboursOfV[neighOfVOffset + bestMove] = 1;
    neighboursOfV[neighOfVOffset + part] = 0;

    for (i = vOffsets[vertex]; i < vertOffset; ++i) {
      // ###
      // update the hyperedge stats: (vInPart etc.)
      // ###

      hEdge = vToHedges[i];
      hEdgeOffset = hEdgeVinPartOffsets[hEdge];

      --hEdgeVinPart[hEdgeOffset + part];
      ++hEdgeVinPart[hEdgeOffset + bestMove];

      if (hEdgeVinPart[hEdgeOffset + part] > 0)
        neighboursOfV[neighOfVOffset + part] = 1;
    }

    if (neighboursOfV[neighOfVOffset + part] == 0)
      --numNeighParts[vertex];

    // ###
    // update the adj vertices stats:
    // (num neighbours in part etc.)
    // ###

    updateAdjVertStat(vertex, part, bestMove);

    // ###
    // update other structs
    // ###

    partitionVector[vertex] = bestMove;
    partWeights[part] -= vertexWt;
    partWeights[bestMove] += vertexWt;

    if (partWeights[part] <= maxPartWt)
      overWeight[part] = 0;

    // ###
    // finally, update the gain...
    // ###

    gain += vGain;
    part = findHeaviestOverweight();
  }

  return gain;
}

#endif
