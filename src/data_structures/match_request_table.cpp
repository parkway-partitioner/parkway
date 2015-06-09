#include "data_structures/match_request_table.hpp"

namespace parkway {
namespace data_structures {


match_request_table::match_request_table(int _size)
    : numEntries(0),
      size(_size) {

  table.setLength(size);
  for (std::size_t i = 0; i < size; ++i) {
    table[i] = nullptr;
  }
}

match_request_table::~match_request_table() {
  for (std::size_t i = 0; i < size; ++i) {
    DynaMem<entry>::deletePtr(table[i]);
  }
}

match_request_table::entry *match_request_table::getEntryPtr(int _vertex) const {
  entry *entry_ = table[_vertex % size];

  while (entry_ && entry_->getNonLocal() != _vertex) {
    entry_ = entry_->getNextEntry();
  }

  return entry_;
}

int match_request_table::lookupClusterWt(int _vertex) const {
  entry *entry_ = table[_vertex % size];

  while (entry_ && entry_->getNonLocal() != _vertex) {
    entry_ = entry_->getNextEntry();
  }

  return entry_ ? entry_->getClusterWt() : -1;
}

int match_request_table::lookupCluIndex(int _vertex) const {
  entry *entry_ = table[_vertex % size];

  while (entry_ && entry_->getNonLocal() != _vertex) {
    entry_ = entry_->getNextEntry();
  }

  return entry_ ? entry_->getCluIndex() : -1;
}

int match_request_table::lookupNumLocals(int _vertex) const {
  entry *entry_ = table[_vertex % size];

  while (entry_ && entry_->getNonLocal() != _vertex) {
    entry_ = entry_->getNextEntry();
  }

  return entry_ ? entry_->getNumLocals() : -1;
}

void match_request_table::clearTable() {
  for (std::size_t i = 0; i < size; ++i)
    DynaMem<entry>::deletePtr(table[i]);

  numEntries = 0;
  entryPtrs.setLength(0);
}

void match_request_table::addLocal(int _vertex, int _local, int locWt,
                                   int proc) {
  int slot = _vertex % size;
  entry *entry_ = table[slot];

  while (entry_ && entry_->getNonLocal() != _vertex) {
    entry_ = entry_->getNextEntry();
  }

  if (entry_) {
    entry_->addLocal(_local, locWt);
  } else {
    entry *newEntry = new entry(_vertex, _local, locWt, proc, table[slot]);
    table[slot] = newEntry;
    entryPtrs.assign(numEntries++, newEntry);
  }
}

void match_request_table::setCluIndex(int vertex, int index, int cluWt) {
  entry *entry_ = table[vertex % size];

  while (entry_ && entry_->getNonLocal() != vertex) {
    entry_ = entry_->getNextEntry();
  }

  if (entry_) {
    entry_->setCluIndex(index);
    entry_->setCluWeight(cluWt);
  }
}

void match_request_table::removeLocal(int vertex, int locVertex, int locWt) {
  entry *entry_ = table[vertex % size];

  while (entry_ && entry_->getNonLocal() != vertex) {
    entry_ = entry_->getNextEntry();
  }

  if (entry_) {
    entry_->removeLocal(locVertex, locWt);
  }
}

void match_request_table::removeEntry(int _vertex) {
  entry *entry_ = table[_vertex % size];

  while (entry_ && entry_->getNonLocal() != _vertex) {
    entry_ = entry_->getNextEntry();
  }

  if (entry_) {
    entry_->clearEntry();
  }
}

}  // namespace data_structures
}  // namespace parkway
