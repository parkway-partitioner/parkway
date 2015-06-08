
#ifndef _BUCKET_HPP
#define _BUCKET_HPP

// ### BucketNode.hpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 30/11/2004: Last Modified
//
// ###

#include "DynamicArray.h"

typedef struct Bucket {
  int vertexID;
  Bucket *prev;
  Bucket *next;

} BucketNode;

typedef BucketNode *BucketNodePtr;
typedef DynamicArray<BucketNode> NodeArray;
typedef DynamicArray<BucketNodePtr> NodePtrArray;

#endif
