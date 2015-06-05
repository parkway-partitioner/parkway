
#  ifndef _HASH_TABLES_HPP
#  define _HASH_TABLES_HPP


// ### HashTables.hpp ###
//
// Copyright (C) 2004 Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
// 
// 28/2/2005: Last Modified
//
// NOTE:
//
// - uses chaining to resolve collisions
//
// ###


#  include "Macros.h"
#  include "Funct.hpp"
#  include "Dyna.hpp"
#  include "CompleteBinaryTree.hpp"


#   define PRIMARY_HASH(a,b) (Mod(a,b))
#   define SECONDARY_HASH(a,b) (1+Mod(a,(b-1)))


using namespace std;


class TableUtils
{
protected:

  static CompleteBinaryTree<int> tableSizeTree;
  static FastDynaArray<int> scatterArray; 
  static int scatterSize;

public:

  TableUtils();
  ~TableUtils() { }
 
  static void setScatterArray(int size);
  
  static inline int scatterKey(register int i) { return scatterArray[i]; }
  static inline int getTableSize(int n) { return (tableSizeTree.getRootVal(n)); }

 
};




/* New MapFromPosInt Class */
								
template <class T>
class MapFromPosInt
{
  /* requires positive int keys */

protected:

  int numEntries;
  int size;

  FastDynaArray<T> table;
  FastDynaArray<int> keys;  

public:

  MapFromPosInt();
  MapFromPosInt(int _size);
  ~MapFromPosInt() { }

  void createTable(int _size);
  void destroyTable();
  void recoverTable();

  int insertKey(int key, T val);
  T &getVal(int key);

  inline int getNumEntries() { return numEntries; }
  inline int getNumSlots() { return size; }
};


/* map from +ve int to +ve int */


class MapToPosInt
{
  /* requires positive int keys */

protected:

  int numEntries;
  int size;
  int useHash;

  FastDynaArray<int> entries;
  FastDynaArray<int> table;
  FastDynaArray<int> keys;  
  
public:

  MapToPosInt();
  MapToPosInt(int _size, int use_hash);
  ~MapToPosInt() { }

  void createTable(int _size, int use_hash);
  void destroyTable();
  void recoverTable();

  int insertKey(int key, int val);
  //int resetKey(int key);
  void resetSlots();
  int insertIfEmpty(int key, int val);
  
  int getCareful(int key);
  int getVal(int key);

  inline int usingHashTable() { return useHash; }
  inline int getNumEntries() { return numEntries; }
  inline int getNumSlots() { return size; }
};




/* new HedgeIndexTable */


class NewHedgeIndexTable
{
protected:

  int numEntries;
  int size;

  FastDynaArray<int> table;
  FastDynaArray<int> nextSameKey;
  FastDynaArray<HashKey> keys;  

public:

  NewHedgeIndexTable(int _size);
  ~NewHedgeIndexTable() { }

  void createTable(int _size);
  void destroyTable();
  void recoverTable();

  void print();
  void insertKey(HashKey key, int index);
  int getHedgeIndex(HashKey key, int &numSeen);

  inline int getNumEntries() { return numEntries; }
  inline int getNumSlots() { return size; }
  
};




//////////////////////////////////////////////
// class declarations for MatchRequestTable //
//////////////////////////////////////////////


class MatchRequestEntry 
{
  
protected:

#  ifdef DEBUG_TABLES
  static int numInitialised;
  static int numDestroyed;
#  endif  
  
  int nonLocVertex;
  int clusterWeight;
  int clusterIndex;
  int nonLocProc;
  int numLocals;

  FastDynaArray<int> locVertices;
  MatchRequestEntry *next;

public:

  inline MatchRequestEntry(int _nonlocal, int _local, int _locWt, int _proc, MatchRequestEntry* _next)
  {
    nonLocVertex = _nonlocal;
    clusterWeight = _locWt;
    clusterIndex = -1;
    nonLocProc = _proc;
    locVertices.assign(0,_local);
    numLocals = 1;
    next = _next;

#  ifdef DEBUG_TABLES
    ++numInitialised;
#  endif
  }
  
  inline ~MatchRequestEntry() 
  {
    DynaMem<MatchRequestEntry>::deletePtr(next);

#  ifdef DEBUG_TABLES
    ++numDestroyed;
#  endif
  }
  
  inline int getNonLocal() const { return nonLocVertex; }
  inline int getClusterWt() const { return clusterWeight; }
  inline int getCluIndex() const { return clusterIndex; }
  inline int getNonLocProc() const { return nonLocProc; }
  inline int getNumLocals() const { return numLocals; }
  inline int *getLocalsArray() const { return locVertices.getArray(); }
  inline MatchRequestEntry *getNextEntry() const { return next; }

  inline void setNextEntry(register MatchRequestEntry* newNext) { next = newNext; }
  inline void setCluIndex(register int _index) { clusterIndex = _index; }
  inline void setNonLocProc(register int _proc) { nonLocProc = _proc; }
  inline void setCluWeight(register int _cluWt) { clusterWeight = _cluWt; }
  inline void clearEntry() { numLocals=0; clusterWeight=0; } 

  inline void addLocal(register int _loc, register int _locWt) 
  {
    locVertices.assign(numLocals++, _loc);
    clusterWeight += _locWt;
  }

  inline void removeLocal(register int loc, int locWt) 
  {
    register int i=0;
    register int j = numLocals-1;

    for ( ;i<numLocals;++i)
      if(locVertices[i] == loc)
	break;

#  ifdef DEBUG_TABLES
    assert(i < numLocals);
#  endif

    while (i < j) 
      { 
	locVertices[i] = locVertices[i+1];
	++i;
      }
    
    --numLocals;
    clusterWeight -= locWt;
  }

#  ifdef DEBUG_TABLES  
  static inline int getNumDeleted() { return numDestroyed; }
  static inline int getNumAllocated() { return numInitialised; }  
#  endif

};


class MatchRequestTable 
{

protected:

  int numEntries;
  int size;
 
  FastDynaArray<MatchRequestEntry*> table;
  FastDynaArray<MatchRequestEntry*> entryPtrs;

public:

  MatchRequestTable(int _size);
  ~MatchRequestTable();

  int lookupClusterWt(register int _vertex) const;
  int lookupCluIndex(register int _vertex) const;
  int lookupNumLocals(register int _vertex) const;
  
  inline int getNumEntries() const { return numEntries; }
  inline int getNumSlots() const { return size; }
  inline MatchRequestEntry *getEntry(register int i) const { return entryPtrs[i]; }
  inline MatchRequestEntry **getEntriesArray() const { return entryPtrs.getArray(); }

  MatchRequestEntry *getEntryPtr(register int _vertex) const; 

  void addLocal(int _vertex, int _local, int locWt, int nonLocProc);
  void setCluIndex(register int _vertex, int _index, int _cluWt);
  void removeLocal(register int _vertex, int _local, int locWt);
  void removeEntry(register int _vertex);
  void clearTable();
};



#  endif
