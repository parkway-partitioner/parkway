#ifndef _READER_CPP
#define _READER_CPP

// ### reader.cpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 3/12/2004: Last Modified
//
// ###

#include "reader.h"

using parkway::hypergraph::parallel_hypergraph;

void initGraphStructs(int &numLocalVertices, int &numLocalHedges,
                      int *&vWeights, int *&hEdgeWts, int *&pinList,
                      int *&offsets, const char *filename, int myRank) {
  int buffer[3];
  int hEdgeDataLen;

  int i;
  int j;
  int ij;

  char my_file[512];
  char message[512];

  ifstream in_stream;

  dynamic_array<int> data;

  sprintf(my_file, "%s-%d", filename, myRank);

  in_stream.open(my_file, ifstream::in | ifstream::binary);

  if (!in_stream.is_open()) {
    sprintf(message, "p[%d] could not open file %s\n", myRank, my_file);
    cout << message;
    MPI_Abort(MPI_COMM_WORLD, 0);
  }

  i = sizeof(int) * 3;
  in_stream.read((char *)(&buffer[0]), i);
  if (in_stream.gcount() != i) {
    sprintf(message, "p[%d] cannot read int buffer[3]\n", myRank);
    cout << message;
    MPI_Abort(MPI_COMM_WORLD, 0);
  }

  numLocalVertices = buffer[1];
  hEdgeDataLen = buffer[2];

  vWeights = new int[numLocalVertices];
  if (!vWeights)
    return (error(myRank, numLocalVertices, "[loc vertices]"));
  data.reserve(hEdgeDataLen);

  i = sizeof(int) * numLocalVertices;
  in_stream.read((char *)(vWeights), i);
  if (in_stream.gcount() != i) {
    sprintf(message, "p[%d] cannot read %d vertices\n", myRank,
            numLocalVertices);
    cout << message;
    MPI_Abort(MPI_COMM_WORLD, 0);
  }

  i = sizeof(int) * hEdgeDataLen;
  in_stream.read((char *)(data.data()), i);
  if (in_stream.gcount() != i) {
    sprintf(message, "p[%d] cannot read %d of hyperedge data_\n", myRank,
            hEdgeDataLen);
    cout << message;
    MPI_Abort(MPI_COMM_WORLD, 0);
  }

  in_stream.close();

  ij = 0;
  for (i = 0; i < hEdgeDataLen; i += data[i])
    ++ij;

  numLocalHedges = ij;
  hEdgeWts = new int[ij];
  if (!hEdgeWts)
    error(myRank, ij, "[h weights]");
  offsets = new int[ij + 1];
  if (!offsets)
    error(myRank, ij + 1, "[h offsets]");

  j = 0;
  offsets[j] = 0;

  for (i = 0; i < hEdgeDataLen; i += ij) {
    ij = data[i];
    hEdgeWts[j] = data[i + 1];
    auto offset = offsets[j++] + (ij - 2);
    offsets[j + 1] = offset;
  }

  pinList = new int[offsets[j]];
  if (!pinList)
    error(myRank, offsets[j], "[pin list]");

  ij = 0;
  for (i = 0; i < hEdgeDataLen; i += data[i])
    for (j = 2; j < data[i]; ++j)
      pinList[ij++] = data[i + j];
}

void error(int myRank, int a, const char *note) {
  char message[512];

  sprintf(message, "p[%d] (%s) memory allocation of size %d failed\n", myRank,
          note, a);
  cout << message;
  MPI_Abort(MPI_COMM_WORLD, 0);
}

void testRecordedPartition(const char *filename, int myRank, int numProcs,
                           int numParts, double constraint, ostream &out,
                           MPI_Comm comm) {
  parallel_hypergraph *h =
      new parallel_hypergraph(myRank, numProcs, filename, 1, out, comm);

  char pFile[512];
  sprintf(pFile, "%s.part.%d", filename, numParts);

  h->initalize_partition_from_file(pFile, numParts, out, comm);
  h->check_partitions(numParts, constraint, out, comm);

  if (h)
    delete h;
}

void testRecordedPartition(const char *filename, const int *pVector,
                           int numLocVerts, int myRank, int numProcs,
                           int numParts, double constraint, ostream &out,
                           MPI_Comm comm) {
  parallel_hypergraph *h =
      new parallel_hypergraph(myRank, numProcs, filename, 1, out, comm);

  h->set_number_of_partitions(1);
  h->copy_in_partition(pVector, numLocVerts, 0);
  h->check_partitions(numParts, constraint, out, comm);

  if (h)
    delete h;
}

#endif
