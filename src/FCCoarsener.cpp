
#ifndef _FCCOARSENER_CPP
#define _FCCOARSENER_CPP

// ### FCCoarsener.cpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 4/1/2005: Last Modified
//
// ###

#include "FCCoarsener.hpp"

FCCoarsener::FCCoarsener(int _min, int _maxwt, double r, int fanOut, int dbWt,
                         int dL)
    : Coarsener(_min, _maxwt, r, dL) {
  utilFanOut = fanOut;
  divByCluWt = dbWt;
}

FCCoarsener::~FCCoarsener() {}

void FCCoarsener::dispCoarsenerOptions(ostream &out) const {
  switch (dispOption) {
  case SILENT:
    break;

  default:

    out << "|- FCC:"
        << " r = " << reductionRatio << " min = " << minNodes
        << " util = " << utilFanOut << endl
        << "|" << endl;
    break;
  }
}

Hypergraph *FCCoarsener::coarsen(const Hypergraph &h) {
  loadHypergraphForCoarsening(h);

  // ###
  // back off if the hypergraph too small
  // ###

  if (numVertices < minNodes)
    return NULL;

  int numNotMatched = numVertices;
  int bestMatch;
  int bestMatchWt = -1;
  int numNeighbours;
  int limitOnIndex;

  int candVertex;
  int candVertexEntry;

  int i;
  int j;
  int ij;

  int v;
  int coarseIndex;

  int endVerOffset;
  int endHedgeOffset;
  int hEdge;

  double reducedBy;
  double maxMatchMetric;
  double metricVal;

  FastDynaArray<int> neighVerts;
  FastDynaArray<int> neighPairWts;
  FastDynaArray<int> vertices(numVertices);
  FastDynaArray<int> vertexAdjEntry(numVertices);
  FastDynaArray<int> *coarseWts = new FastDynaArray<int>;

  FastDynaArray<double> connectVals;

  BitField toLoad(numHedges);
  toLoad.set1();

  for (i = 0; i < numVertices; ++i) {
    vertices[i] = i;
    vertexAdjEntry[i] = -1;

#ifdef DEBUG_COARSENER
    assert(matchVector[i] == -1);
#endif
  }

  if (currPercentile < 100)
    computeHedgesToLoad(toLoad);

  Funct::randomPermutation(vertices.getArray(), numVertices);

  metricVal = static_cast<double>(numVertices) / reductionRatio;
  limitOnIndex = numVertices - static_cast<int>(floor(metricVal - 1.0));

  i = 0;
  coarseIndex = 0;

  for (; i < numVertices; ++i) {
    if (matchVector[vertices[i]] == -1) {
      v = vertices[i];
      bestMatch = -1;
      maxMatchMetric = 0;
      numNeighbours = 0;

      endVerOffset = vOffsets[v + 1];

      for (j = vOffsets[v]; j < endVerOffset; ++j) {
        hEdge = vToHedges[j];

        if (toLoad(hEdge) == 1) {
          endHedgeOffset = hEdgeOffsets[hEdge + 1];

          for (ij = hEdgeOffsets[hEdge]; ij < endHedgeOffset; ++ij) {
            candVertex = pinList[ij];

            if (candVertex != v) {
              candVertexEntry = vertexAdjEntry[candVertex];

              if (candVertexEntry == -1) {
                neighVerts.assign(numNeighbours, candVertex);

                if (matchVector[candVertex] == -1)
                  neighPairWts.assign(numNeighbours,
                                      vWeight[v] + vWeight[candVertex]);
                else
                  neighPairWts.assign(
                      numNeighbours,
                      vWeight[v] + (*coarseWts)[matchVector[candVertex]]);

                if (utilFanOut)
                  connectVals.assign(
                      numNeighbours,
                      static_cast<double>(hEdgeWeight[hEdge]) /
                          (endHedgeOffset - (hEdgeOffsets[hEdge] + 1)));
                else
                  connectVals.assign(numNeighbours,
                                     static_cast<double>(hEdgeWeight[hEdge]));

                vertexAdjEntry[candVertex] = numNeighbours;
                ++numNeighbours;
              } else {
                if (utilFanOut)
                  connectVals[candVertexEntry] +=
                      static_cast<double>(hEdgeWeight[hEdge]) /
                      (endHedgeOffset - (hEdgeOffsets[hEdge] + 1));
                else
                  connectVals[candVertexEntry] +=
                      static_cast<double>(hEdgeWeight[hEdge]);
              }
            }
          }
        }
      }

      if (numNeighbours == 0) {
        // match vertex v as singleton as not connected

        matchVector[v] = coarseIndex;
        coarseWts->assign(coarseIndex++, vWeight[v]);
        --numNotMatched;
      } else {
        // else need to pick the best match

        for (j = 0; j < numNeighbours; ++j) {
          if (neighPairWts[j] < maxVertexWt) {
            metricVal = connectVals[j];

            if (divByCluWt)
              metricVal /= neighPairWts[j];

            if (metricVal > maxMatchMetric) {
              maxMatchMetric = metricVal;
              bestMatch = neighVerts[j];
              bestMatchWt = neighPairWts[j];
            }
          }

          vertexAdjEntry[neighVerts[j]] = -1;
        }

        if (bestMatch == -1) {
          // no best match as all would violate max
          // vertex weight - match as singleton

          matchVector[v] = coarseIndex;
          coarseWts->assign(coarseIndex++, vWeight[v]);
          --numNotMatched;
        } else {
          if (matchVector[bestMatch] == -1) {
            // match with another unmatched

            matchVector[v] = coarseIndex;
            matchVector[bestMatch] = coarseIndex;
            coarseWts->assign(coarseIndex++, bestMatchWt);
            numNotMatched -= 2;
          } else {
            // match with existing cluster of coarse vertices

            matchVector[v] = matchVector[bestMatch];
            (*coarseWts)[matchVector[v]] += vWeight[v];
            --numNotMatched;
          }
        }
      }

      // check if hypergraph sufficiently shrunk

      if (i > limitOnIndex) {
        reducedBy =
            static_cast<double>(numVertices) / (numNotMatched + coarseIndex);

        if (reducedBy > reductionRatio)
          break;
      }
    }
  }

  if (i == numVertices) {
    // cannot sufficiently reduce hypergraph
    // back off

    DynaMem<FastDynaArray<int> >::deletePtr(coarseWts);
    return NULL;
  }

  // match remaining unmatched vertices as singletons

  for (; i < numVertices; ++i) {
    v = vertices[i];

    if (matchVector[v] == -1) {
      matchVector[v] = coarseIndex;
      coarseWts->assign(coarseIndex++, vWeight[v]);
    }
  }

#ifdef DEBUG_COARSENER
  for (j = 0; j < numVertices; ++j)
    assert(matchVector[j] >= 0 && matchVector[j] < coarseIndex);
#endif

  coarseWts->setLength(coarseIndex);

  return (buildCoarseHypergraph(coarseWts->getArray(), coarseIndex,
                                h.getTotWeight()));
}

#endif
