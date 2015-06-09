
#ifndef _BUCKET_HPP
#define _BUCKET_HPP

// ### bucket_node.hpp ###
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

typedef struct bucket {
  int vertexID;
  bucket *prev;
  bucket *next;
} bucket_node;

typedef bucket_node *BucketNodePtr;
typedef dynamic_array<bucket_node> NodeArray;
typedef dynamic_array<BucketNodePtr> NodePtrArray;

}  // namespace data_structures
}  // namespace parkway

#endif
