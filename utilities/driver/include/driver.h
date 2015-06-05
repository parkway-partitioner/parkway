

#ifndef _DRIVER_H
#define _DRIVER_H

// ### driver.h ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 26/11/2004: Created
//
// ###

#include <fstream>
#include "mpi.h"
#include "reader.h"

using namespace std;

void ParaPartKway(const char *file_name, const char *outFile, int num_parts,
                  double constraint, int &k_1cut, const int *options,
                  MPI_Comm comm);
void ParaPartKway(int numVertices, int numHedges, const int *vWeights,
                  const int *hEdgeWts, const int *pinList, const int *offsets,
                  int numParts, double constraint, int &k_1cut,
                  const int *options, int *pVector, const char *outFile,
                  MPI_Comm comm);

#endif
