
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

#include "data_structures/dynamic_array.hpp"

namespace parkway {
namespace data_structures {

typedef struct Bucket {
  int vertexID;
  Bucket *prev;
  Bucket *next;
} BucketNode;

typedef BucketNode *BucketNodePtr;
typedef dynamic_array<BucketNode> NodeArray;
typedef dynamic_array<BucketNodePtr> NodePtrArray;

}  // namespace data_structures
}  // namespace parkway

#endif
