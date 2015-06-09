#ifndef _DATA_STRUCTURES_MATCH_REQUEST_TABLE_HPP
#define _DATA_STRUCTURES_MATCH_REQUEST_TABLE_HPP

#include "data_structures/dynamic_array.hpp"

namespace parkway {
namespace data_structures {

class match_request_table {
 public:
  class entry;

  match_request_table(int _size);
  ~match_request_table();

  int lookupClusterWt(int _vertex) const;
  int lookupCluIndex(int _vertex) const;
  int lookupNumLocals(int _vertex) const;

  inline int getNumEntries() const { return numEntries; }
  inline int getNumSlots() const { return size; }
  inline entry *getEntry(int i) const {
    return entryPtrs[i];
  }
  inline entry **getEntriesArray() const {
    return entryPtrs.getArray();
  }

  entry *getEntryPtr(int _vertex) const;

  void addLocal(int _vertex, int _local, int locWt, int nonLocProc);
  void setCluIndex(int _vertex, int _index, int _cluWt);
  void removeLocal(int _vertex, int _local, int locWt);
  void removeEntry(int _vertex);
  void clearTable();

 protected:
  int numEntries;
  int size;

  DynamicArray<entry *> table;
  DynamicArray<entry *> entryPtrs;
};


class match_request_table::entry {
 protected:
  int nonLocVertex;
  int clusterWeight;
  int clusterIndex;
  int nonLocProc;
  int numLocals;

  DynamicArray<int> locVertices;
  entry *next;

 public:
  inline entry(int _nonlocal, int _local, int _locWt, int _proc, entry *_next)
      : nonLocVertex(_nonlocal),
        clusterWeight(_locWt),
        clusterIndex(-1),
        nonLocProc(_proc),
        numLocals(1),
        next(_next) {
    locVertices.assign(0, _local);
  }

  inline ~entry() {
    DynaMem<entry>::deletePtr(next);
  }

  inline int getNonLocal() const { return nonLocVertex; }
  inline int getClusterWt() const { return clusterWeight; }
  inline int getCluIndex() const { return clusterIndex; }
  inline int getNonLocProc() const { return nonLocProc; }
  inline int getNumLocals() const { return numLocals; }
  inline int *getLocalsArray() const { return locVertices.getArray(); }
  inline entry *getNextEntry() const { return next; }

  inline void setNextEntry(entry *newNext) {
    next = newNext;
  }
  inline void setCluIndex(int _index) { clusterIndex = _index; }
  inline void setNonLocProc(int _proc) { nonLocProc = _proc; }
  inline void setCluWeight(int _cluWt) { clusterWeight = _cluWt; }
  inline void clearEntry() {
    numLocals = 0;
    clusterWeight = 0;
  }

  inline void addLocal(int _loc, int _locWt) {
    locVertices.assign(numLocals++, _loc);
    clusterWeight += _locWt;
  }

  inline void removeLocal(int loc, int locWt) {
    int i = 0;
    int j = numLocals - 1;

    for (; i < numLocals; ++i)
      if (locVertices[i] == loc)
        break;

    while (i < j) {
      locVertices[i] = locVertices[i + 1];
      ++i;
    }

    --numLocals;
    clusterWeight -= locWt;
  }
};

}  // namespace data_structures
}  // namespace parkway

#endif  // _DATA_STRUCTURES_MATCH_REQUEST_TABLE_HPP
