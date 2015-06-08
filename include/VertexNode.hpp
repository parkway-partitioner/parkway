#ifndef _VERTEX_NODE_HPP
#define _VERTEX_NODE_HPP

// ### VertexNode.hpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 30/11/2004: Last Modified
//
// ###

#include "DynamicArray.h"

typedef struct VertexNode {
  int vertexID;
  VertexNode *next;

} VNode;

typedef VNode *VNodePtr;
typedef DynamicArray<VNodePtr> VNodePtrArray;

#endif
