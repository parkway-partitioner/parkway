
#ifndef _HASH_TABLES_CPP
#define _HASH_TABLES_CPP

// ### HashTables.cpp ###
//
// Copyright (C) 2004, Aleksandar Trifunovic, Imperial College London
//
// HISTORY:
//
// 13/4/2004: Added DuplRemEntry and DuplRemTable
//            classes.
//
// 3/2/2005: Last Modified
//
// ###

#include "HashTables.hpp"

CompleteBinaryTree<int> TableUtils::tableSizeTree;
FastDynaArray<int> TableUtils::scatterArray;
int TableUtils::scatterSize = 0;

TableUtils::TableUtils() {
  int sizes[16];
  int mins[16];

  sizes[0] = 1091;
  sizes[1] = 2113;
  sizes[2] = 4133;
  sizes[3] = 9067;
  sizes[4] = 17093;
  sizes[5] = 37097;
  sizes[6] = 70099;
  sizes[7] = 145109;
  sizes[8] = 300149;
  sizes[9] = 610217;
  sizes[10] = 1290151;
  sizes[11] = 2600177;
  sizes[12] = 6000109;
  sizes[13] = 12500197;
  sizes[14] = 26000111;
  sizes[15] = 64000147;

  mins[0] = 0;
  mins[1] = 256;
  mins[2] = 512;
  mins[3] = 1024;
  mins[4] = 2056;
  mins[5] = 4112;
  mins[6] = 8224;
  mins[7] = 16448;
  mins[8] = 35000;
  mins[9] = 70000;
  mins[10] = 140000;
  mins[11] = 300000;
  mins[12] = 700000;
  mins[13] = 1800000;
  mins[14] = 4000000;
  mins[15] = 12500000;

  tableSizeTree.setupTree(sizes, mins, 16);
}

void TableUtils::setScatterArray(int size) {
  scatterSize = size;
  scatterArray.setLength(scatterSize);

  register int i;

  for (i = 0; i < scatterSize; ++i)
    scatterArray[i] = i;

  Funct::randomPermutation(scatterArray.getArray(), scatterSize);
}

/* New Map Class */

template <class T> MapFromPosInt<T>::MapFromPosInt() { createTable(0); }

template <class T> MapFromPosInt<T>::MapFromPosInt(int _size) {
  createTable(_size);
}

template <class T> void MapFromPosInt<T>::createTable(int _size) {
  numEntries = 0;
  size = TableUtils::getTableSize(_size);

#ifdef DEBUG_TABLES
  assert(size >= _size);
#endif

  keys.setLength(size);
  table.setLength(size);

  register int i;

  for (i = 0; i < size; ++i)
    keys[i] = -1;
}

template <class T> void MapFromPosInt<T>::destroyTable() {
  numEntries = 0;
  table.setLength(0);
  keys.setLength(0);
}

template <class T> void MapFromPosInt<T>::recoverTable() {
  table.setLength(size);
  keys.setLength(size);
  numEntries = 0;

  register int i;

  for (i = 0; i < size; ++i)
    keys[i] = -1;
}

template <class T> int MapFromPosInt<T>::insertKey(int key, T val) {
  register int indepKey = TableUtils::scatterKey(key);
  register int slot = PRIMARY_HASH(indepKey, size);

  while (keys[slot] != -1 && keys[slot] != indepKey)
    slot = Mod(slot + SECONDARY_HASH(indepKey, size), size);

  if (keys[slot] == -1) {
    table[slot] = val;
    keys[slot] = indepKey;
    ++numEntries;
    return 0;
  } else {
    table[slot] = val;
    return 1;
  }
}

template <class T> T &MapFromPosInt<T>::getVal(int key) {
  register int indepKey = TableUtils::scatterKey(key);
  register int slot = PRIMARY_HASH(indepKey, size);

  while (keys[slot] != indepKey) {
    slot = Mod(slot + SECONDARY_HASH(indepKey, size), size);
#ifdef DEBUG_TABLES
    assert(keys[slot] != -1);
#endif
  }

  return table[slot];
}

template class MapFromPosInt<int>;

/* map from +ve int to +ve int */

MapToPosInt::MapToPosInt() { createTable(0, 0); }

MapToPosInt::MapToPosInt(int _size, int use_hash) {
  createTable(_size, use_hash);
}

void MapToPosInt::createTable(int _size, int use_hash) {
  useHash = use_hash;
  numEntries = 0;
  entries.setLength(2048);

  register int i;

  if (useHash == 0) {
    size = _size;
    keys.setLength(0);
    table.setLength(size);

    for (i = 0; i < size; ++i)
      table[i] = -1;
  } else {
    size = TableUtils::getTableSize(_size);
#ifdef DEBUG_TABLES
    assert(size >= _size);
#endif

    keys.setLength(size);
    table.setLength(size);

    for (i = 0; i < size; ++i) {
      keys[i] = -1;
      table[i] = -1;
    }
  }
}

void MapToPosInt::destroyTable() {
  numEntries = 0;
  entries.setLength(0);
  table.setLength(0);
  keys.setLength(0);
}

void MapToPosInt::recoverTable() {
  register int i;

  numEntries = 0;
  entries.setLength(2048);

  if (useHash) {
    table.setLength(size);
    keys.setLength(size);

    for (i = 0; i < size; ++i) {
      keys[i] = -1;
      table[i] = -1;
    }
  } else {
    table.setLength(size);
    keys.setLength(0);

    for (i = 0; i < size; ++i)
      table[i] = -1;
  }
}

int MapToPosInt::insertKey(int key, int val) {
#ifdef DEBUG_TABLES
  assert(key >= 0);
  assert(val >= 0);
#endif

  if (useHash) {
    register int indepKey = TableUtils::scatterKey(key);
    register int slot = PRIMARY_HASH(indepKey, size);

    while (keys[slot] != -1 && keys[slot] != indepKey)
      slot = Mod(slot + SECONDARY_HASH(indepKey, size), size);

    if (keys[slot] == -1) {
#ifdef DEBUG_TABLES
      assert(table[slot] == -1);
#endif
      table[slot] = val;
      keys[slot] = indepKey;
      entries.assign(numEntries++, slot);
      return 0;
    } else {
#ifdef DEBUG_TABLES
      assert(keys[slot] == indepKey);
      assert(table[slot] >= 0);
#endif
      table[slot] = val;
      return 1;
    }
  } else {
#ifdef DEBUG_TABLES
    assert(key < size);
#endif
    if (table[key] == -1) {
      table[key] = val;
      entries.assign(numEntries++, key);
      return 0;
    } else {
      table[key] = val;
      return 1;
    }
  }
}

int MapToPosInt::insertIfEmpty(int key, int val) {
#ifdef DEBUG_TABLES
  assert(key >= 0);
#endif

  if (useHash) {
    register int indepKey = TableUtils::scatterKey(key);
    register int slot = PRIMARY_HASH(indepKey, size);

    while (keys[slot] != -1 && keys[slot] != indepKey)
      slot = Mod(slot + SECONDARY_HASH(indepKey, size), size);

    if (keys[slot] == -1) {
      table[slot] = val;
      keys[slot] = indepKey;
      entries.assign(numEntries++, slot);
      return -1;
    } else {
      return table[slot];
    }
  } else {
#ifdef DEBUG_TABLES
    assert(key < size);
#endif

    if (table[key] == -1) {
      table[key] = val;
      entries.assign(numEntries++, key);
      return -1;
    } else {
      return table[key];
    }
  }
}

void MapToPosInt::resetSlots() {
  register int i;
  register int slot;

  if (useHash) {
    for (i = 0; i < numEntries; ++i) {
      slot = entries[i];

      keys[slot] = -1;
      table[slot] = -1;
    }
  } else {
    for (i = 0; i < numEntries; ++i)
      table[entries[i]] = -1;
  }

  numEntries = 0;
}

int MapToPosInt::getCareful(int key) {
#ifdef DEBUG_TABLES
  assert(key >= 0);
#endif

  if (useHash) {
    register int indepKey = TableUtils::scatterKey(key);
    register int slot = PRIMARY_HASH(indepKey, size);

    while (keys[slot] != indepKey && keys[slot] != -1)
      slot = Mod(slot + SECONDARY_HASH(indepKey, size), size);

    if (keys[slot] == -1) {
#ifdef DEBUG_TABLES
      assert(table[slot] == -1);
#endif
      return -1;
    } else {
#ifdef DEBUG_TABLES
      assert(table[slot] != -1);
#endif
      return table[slot];
    }
  } else {
#ifdef DEBUG_TABLES
    assert(key < size);
#endif
    if (table[key] == -1)
      return -1;
    else
      return table[key];
  }
}

int MapToPosInt::getVal(int key) {
#ifdef DEBUG_TABLES
  assert(key >= 0);
#endif

  if (useHash) {
    register int indepKey = TableUtils::scatterKey(key);
    register int slot = PRIMARY_HASH(indepKey, size);

    while (keys[slot] != indepKey) {
      slot = Mod(slot + SECONDARY_HASH(indepKey, size), size);
#ifdef DEBUG_TABLES
      assert(keys[slot] != -1);
#endif
    }

    return table[slot];
  } else {
#ifdef DEBUG_TABLES
    assert(key < size);
#endif
    return table[key];
  }
}

/* class NewHedgeIndexTable */

NewHedgeIndexTable::NewHedgeIndexTable(int _size) { createTable(_size); }

void NewHedgeIndexTable::createTable(int _size) {
  numEntries = 0;
  size = TableUtils::getTableSize(_size);

#ifdef DEBUG_TABLES
  assert(size >= _size);
#endif

  keys.setLength(size);
  table.setLength(size);
  nextSameKey.setLength(size);

  register int i;

  for (i = 0; i < size; ++i) {
    table[i] = -1;
    nextSameKey[i] = -1;
  }
}

void NewHedgeIndexTable::destroyTable() {
  numEntries = 0;
  table.setLength(0);
  keys.setLength(0);
  nextSameKey.setLength(0);
}

void NewHedgeIndexTable::recoverTable() {
  table.setLength(size);
  keys.setLength(size);
  nextSameKey.setLength(size);
  numEntries = 0;

  register int i;

  for (i = 0; i < size; ++i) {
    table[i] = -1;
    nextSameKey[i] = -1;
  }
}

void NewHedgeIndexTable::insertKey(HashKey key, int index) {
  register int slot = PRIMARY_HASH(key, size);
  register int lastSeen = -1;

  while (table[slot] != -1) {
    if (keys[slot] == key)
      lastSeen = slot;

    slot = Mod(slot + SECONDARY_HASH(key, size), size);
  }

  table[slot] = index;
  keys[slot] = key;

  if (lastSeen >= 0)
    nextSameKey[lastSeen] = slot;

  ++numEntries;
#ifdef DEBUG_TABLES
  assert(numEntries < size);
#endif
}

int NewHedgeIndexTable::getHedgeIndex(HashKey key, int &numSeen) {
#ifdef DEBUG_TABLES
  assert(numSeen < size);
#endif

  register int slot;

  if (numSeen == -1) {
    slot = PRIMARY_HASH(key, size);

    while (table[slot] != -1) {
      if (keys[slot] == key) {
        numSeen = nextSameKey[slot];
#ifdef DEBUG_TABLES
        assert(numSeen < size);
#endif
        return table[slot];
      }

      slot = Mod(slot + SECONDARY_HASH(key, size), size);
    }

    return -1;
  }

  slot = numSeen;
  numSeen = nextSameKey[slot];
#ifdef DEBUG_TABLES
  assert(numSeen < size);
#endif
  return table[slot];
}

void NewHedgeIndexTable::print() {
  for (int i = 0; i < size; ++i) {
    if (table[i] != -1)
      cout << "[" << i << "]: "
           << "key = " << keys[i] << " index = " << table[i]
           << " next = " << nextSameKey[i] << endl;
    else
      cout << "[empty]" << endl;
  }
}

////////////////////////////////////////////////
// function definitions for MatchRequestTable //
////////////////////////////////////////////////

#ifdef DEBUG_TABLES
int MatchRequestEntry::numInitialised = 0;
int MatchRequestEntry::numDestroyed = 0;
#endif

MatchRequestTable::MatchRequestTable(int _size) {
  size = _size;
  numEntries = 0;

  table.setLength(size);

  register int i = 0;

  for (; i < size; ++i)
    table[i] = NULL;
}

MatchRequestTable::~MatchRequestTable() {
  register int i = 0;

  for (; i < size; ++i)
    DynaMem<MatchRequestEntry>::deletePtr(table[i]);
}

MatchRequestEntry *MatchRequestTable::getEntryPtr(register int _vertex) const {
  register MatchRequestEntry *entry = table[Mod(_vertex, size)];

  while (entry && entry->getNonLocal() != _vertex)
    entry = entry->getNextEntry();

  return entry;
}

int MatchRequestTable::lookupClusterWt(register int _vertex) const {
  register MatchRequestEntry *entry = table[Mod(_vertex, size)];

  while (entry && entry->getNonLocal() != _vertex)
    entry = entry->getNextEntry();

  if (entry)
    return entry->getClusterWt();
  else
    return -1;
}

int MatchRequestTable::lookupCluIndex(register int _vertex) const {
  register MatchRequestEntry *entry = table[Mod(_vertex, size)];

  while (entry && entry->getNonLocal() != _vertex)
    entry = entry->getNextEntry();

  if (entry)
    return entry->getCluIndex();
  else
    return -1;
}

int MatchRequestTable::lookupNumLocals(register int _vertex) const {
  register MatchRequestEntry *entry = table[Mod(_vertex, size)];

  while (entry && entry->getNonLocal() != _vertex)
    entry = entry->getNextEntry();

  if (entry)
    return entry->getNumLocals();
  else
    return -1;
}

void MatchRequestTable::clearTable() {
  register int i = 0;

  for (; i < size; ++i)
    DynaMem<MatchRequestEntry>::deletePtr(table[i]);

  numEntries = 0;
  entryPtrs.setLength(0);
}

void MatchRequestTable::addLocal(int _vertex, int _local, int locWt, int proc) {
  register int slot = Mod(_vertex, size);
  register MatchRequestEntry *entry = table[slot];

  while (entry && entry->getNonLocal() != _vertex)
    entry = entry->getNextEntry();

  if (entry)
    entry->addLocal(_local, locWt);

  else {
    MatchRequestEntry *newEntry =
        new MatchRequestEntry(_vertex, _local, locWt, proc, table[slot]);
    table[slot] = newEntry;
    entryPtrs.assign(numEntries++, newEntry);
  }
}

void MatchRequestTable::setCluIndex(register int vertex, int index, int cluWt) {
  register MatchRequestEntry *entry = table[Mod(vertex, size)];

  while (entry && entry->getNonLocal() != vertex)
    entry = entry->getNextEntry();

#ifdef DEBUG_TABLES
  assert(entry);
#endif

  if (entry) {
    entry->setCluIndex(index);
    entry->setCluWeight(cluWt);
  }
}

void MatchRequestTable::removeLocal(register int vertex, int locVertex,
                                    int locWt) {
  register MatchRequestEntry *entry = table[Mod(vertex, size)];

  while (entry && entry->getNonLocal() != vertex)
    entry = entry->getNextEntry();

#ifdef DEBUG_TABLES
  assert(entry);
#endif

  if (entry)
    entry->removeLocal(locVertex, locWt);
}

void MatchRequestTable::removeEntry(register int _vertex) {
  register MatchRequestEntry *entry = table[Mod(_vertex, size)];

  while (entry && entry->getNonLocal() != _vertex)
    entry = entry->getNextEntry();

  if (entry)
    entry->clearEntry();
}

#endif
