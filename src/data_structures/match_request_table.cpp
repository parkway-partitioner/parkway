#include "data_structures/match_request_table.hpp"

namespace parkway {
namespace data_structures {

match_request_table::match_request_table(int _size)
    : size_(0),
      capacity_(_size) {
  table_.reserve(capacity_);
  for (std::size_t i = 0; i < capacity_; ++i) {
    table_[i] = nullptr;
  }
}

match_request_table::~match_request_table() {
}

match_request_table::entry *match_request_table::get_entry(int _vertex) const {
  entry *entry_ = table_[_vertex % capacity_];
  while (entry_ && entry_->non_local_vertex() != _vertex) {
    entry_ = entry_->next();
  }

  return entry_;
}

int match_request_table::cluster_weight(int _vertex) const {
  entry *entry_ = table_[_vertex % capacity_];

  while (entry_ && entry_->non_local_vertex() != _vertex) {
    entry_ = entry_->next();
  }

  return entry_ ? entry_->cluster_weight() : -1;
}

int match_request_table::cluster_index(int _vertex) const {
  entry *entry_ = table_[_vertex % capacity_];

  while (entry_ && entry_->non_local_vertex() != _vertex) {
    entry_ = entry_->next();
  }

  return entry_ ? entry_->cluster_index() : -1;
}

int match_request_table::local_count(int _vertex) const {
  entry *entry_ = table_[_vertex % capacity_];

  while (entry_ && entry_->non_local_vertex() != _vertex) {
    entry_ = entry_->next();
  }

  return entry_ ? entry_->number_local() : -1;
}

void match_request_table::clear() {
  size_ = 0;
  entries_.reserve(0);
}

void match_request_table::add_local(int _vertex, int _local, int locWt,
                                    int proc) {
  int slot = _vertex % capacity_;
  entry *entry_ = table_[slot];

  while (entry_ && entry_->non_local_vertex() != _vertex) {
    entry_ = entry_->next();
  }

  if (entry_) {
    entry_->add_local(_local, locWt);
  } else {
    entry *newEntry = new entry(_vertex, _local, locWt, proc, table_[slot]);
    table_[slot] = newEntry;
    entries_.assign(size_++, newEntry);
  }
}

void match_request_table::set_cluster_index(int vertex, int index, int cluWt) {
  entry *entry_ = table_[vertex % capacity_];

  while (entry_ && entry_->non_local_vertex() != vertex) {
    entry_ = entry_->next();
  }

  if (entry_) {
    entry_->set_cluster_index(index);
    entry_->set_cluster_weight(cluWt);
  }
}

void match_request_table::remove_local(int vertex, int locVertex, int locWt) {
  entry *entry_ = table_[vertex % capacity_];

  while (entry_ && entry_->non_local_vertex() != vertex) {
    entry_ = entry_->next();
  }

  if (entry_) {
    entry_->remove_local(locVertex, locWt);
  }
}

void match_request_table::remove_entry(int _vertex) {
  entry *entry_ = table_[_vertex % capacity_];

  while (entry_ && entry_->non_local_vertex() != _vertex) {
    entry_ = entry_->next();
  }

  if (entry_) {
    entry_->clear();
  }
}

}  // namespace data_structures
}  // namespace parkway
