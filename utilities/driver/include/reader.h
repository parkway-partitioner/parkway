#ifndef _READER_H
#define _READER_H

// ### reader.h ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 1/12/2004: Last Modified
//
// ###

#include <iostream>
#include <cstdio>
#include "hypergraph/parallel/hypergraph.hpp"

namespace parallel = parkway::parallel;

using namespace std;

void initGraphStructs(int &numLocalVertices, int &numLocalHedges,
                      int *&vWeights, int *&hEdgeWts, int *&pinList,
                      int *&offsets, const char *filename, int myRank);
void error(int myRank, int a, const char *note);

void testRecordedPartition(const char *filename, int myRank, int numProcs,
                           int numParts, double constraint, MPI_Comm comm);
void testRecordedPartition(const char *filename, const int *pVector,
                           int numLocVerts, int myRank, int numProcs,
                           int numParts, double constraint, MPI_Comm comm);

#endif
