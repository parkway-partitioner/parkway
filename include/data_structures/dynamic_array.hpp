#ifndef DATA_STRUCTURES_DYNAMIC_ARRAY_HPP_
#define DATA_STRUCTURES_DYNAMIC_ARRAY_HPP_
#include <vector>
#include <iostream>
#include <iterator>
#include <memory>
#include "utility/sorting.hpp"
#include "utility/random.hpp"

namespace parkway {
namespace data_structures {

template <typename T> class dynamic_array {
 private:
  typedef T value_type;
  typedef std::vector<value_type> data_type;
  typedef typename data_type::size_type size_type;

  // Shared point to Underlying data.
  std::shared_ptr<data_type> data_;

 public:
  // public typedefs.
  typedef value_type& reference;
  typedef const value_type& const_reference;
  typedef value_type&& move_reference;
  typedef typename data_type::iterator iterator;
  typedef typename data_type::const_iterator const_iterator;
  typedef typename data_type::reverse_iterator reverse_iterator;
  typedef typename data_type::const_reverse_iterator const_reverse_iterator;

  // Constructors
  dynamic_array() : data_(std::make_shared<data_type>()) {
  }

  explicit dynamic_array(size_type capacity, const_reference value)
      : data_(std::make_shared<data_type>(capacity, value)) {
  }

  explicit dynamic_array(size_type capacity)
      : data_(std::make_shared<data_type>(capacity)) {
  }

  dynamic_array(dynamic_array& other) {
    data_ = other.data_;
  }

  dynamic_array(const dynamic_array& other) {
    data_ = other.data_;
  }

  dynamic_array& operator=(const dynamic_array& other) {
    data_ = other.data_;
    return *this;
  }

  dynamic_array& operator=(dynamic_array&& other) {
    data_ = other.data_;
    return *this;
  }

  inline void assign(size_type count, const_reference value) {
    data_->assign(count, value);
  }


  // Access
  inline reference at(size_type position) {
    if (this->size() <= position) {
      this->resize(position + 1);
    }
    return data_->at(position);
  }

  inline const_reference at(size_type position) const {
    return data_->at(position);
  }

  inline reference operator[](size_type position) {
    if (this->size() <= position) {
      this->resize(position + 1);
    }
    return data_->operator[](position);
  }

  inline const_reference operator[](size_type position) const {
    return data_->operator[](position);
  }

  inline reference front() {
    return data_->front();
  }

  inline const_reference front() const {
    return data_->front();
  }

  inline reference back() {
    return data_->back();
  }

  inline const_reference back() const {
    return data_->back();
  }

  inline value_type* data() {
    return data_->data();
  }

  inline const value_type* data() const {
    return data_->data();
  }


  // Iterators
  inline iterator begin() {
    return data_->begin();
  }

  inline const_iterator begin() const {
    return data_->begin();
  }

  inline const_iterator cbegin() const {
    return data_->cbegin();
  }

  inline iterator end() {
    return data_->end();
  }

  inline const_iterator end() const {
    return data_->end();
  }

  inline const_iterator cend() const {
    return data_->cend();
  }

  inline reverse_iterator rbegin() {
    return data_->rbegin();
  }

  inline const_reverse_iterator rbegin() const {
    return data_->rbegin();
  }

  inline const_reverse_iterator crbegin() const {
    return data_->crbegin();
  }

  inline reverse_iterator rend() {
    return data_->rend();
  }

  inline const_reverse_iterator rend() const {
    return data_->rend();
  }

  inline const_reverse_iterator crend() const {
    return data_->crend();
  }


  // Capacity
  inline bool empty() const {
    return data_->empty();
  }

  inline size_type size() const {
    return data_->size();
  }

  inline size_type capacity() const {
    return data_->capacity();
  }

  inline void reserve(size_type new_capacity) {
    data_->reserve(new_capacity);
  }

  inline size_type max_size() const {
    return data_->max_size();
  }

  inline void shrink_to_fit() {
    data_->shrink_to_fit();
  }


  // Modifiers
  inline void clear() {
    data_->clear();
  }

  inline void clear_and_shrink() {
    this->clear();
    this->shrink_to_fit();
  }

#if ((__GNUC__ == 4) && (__GNUC_MINOR_ < 9))
  inline iterator insert(iterator position, const_reference value) {
    return data_->insert(position, value);
  }

  inline iterator insert(iterator position, move_reference value) {
    return data_->insert(position, value);
  }

  inline iterator insert(iterator position, size_type count,
                         const_reference value) {
    return data_->insert(position, count, value);
  }

  template <class InputIt>
  inline iterator insert(iterator position, InputIt first, InputIt last) {
    return data_->insert(position, first, last);
  }

  inline iterator insert(iterator position,
                         std::initializer_list<value_type> initializer_list) {
    return data_->insert(position, initializer_list);
  }
#else

  inline iterator insert(const_iterator position, const_reference value) {
    return data_->insert(position, value);
  }

  inline iterator insert(const_iterator position, move_reference value) {
    return data_->insert(position, value);
  }

  inline iterator insert(const_iterator position, size_type count,
                         const_reference value) {
    return data_->insert(position, count, value);
  }

  template <class InputIt>
  inline iterator insert(const_iterator position, InputIt first, InputIt last) {
    return data_->insert(position, first, last);
  }

  inline iterator insert(const_iterator position,
                         std::initializer_list<value_type> initializer_list) {
    return data_->insert(position, initializer_list);
  }
#endif

  inline void set_data(value_type *array, size_type length) {
    this->clear();
    data_->insert(this->begin(), array, array + length);
  }

#if ((__GNUC__ == 4) && (__GNUC_MINOR_ < 9))
  inline iterator erase(iterator position) {
    return data_->erase(position);
  }

  inline iterator erase(iterator first, iterator last) {
    return data_->erase(first, last);
  }
#else
  inline iterator erase(const_iterator position) {
    return data_->erase(position);
  }

  inline iterator erase(const_iterator first, const_iterator last) {
    return data_->erase(first, last);
  }
#endif

  inline void push_back(const_reference value) {
    data_->push_back(value);
  }

  inline void push_back(move_reference value) {
    data_->push_back(value);
  }

  inline void pop_back() {
    data_->pop_back();
  }

  inline void resize(size_type count) {
    data_->resize(count);
  }

  inline void resize(size_type count, const_reference value) {
    data_->resize(count, value);
  }

  inline bool read_from(std::istream &stream, const std::size_t length) {
    if (stream.good()) {
      if (this->size() < length) {
        this->resize(length);
      }
      std::size_t buffer_size = sizeof(value_type) * length;
      stream.read((char *) this->data(), buffer_size);
      return stream.gcount() == buffer_size;
    }
    return false;
  }

  inline void sort() {
    this->sort_between(0, this->size() - 1);
  }

  inline void sort_between(std::size_t start, std::size_t end) {
    parkway::utility::quick_sort(start, end, this->data());
  }

  inline void sort_using_another_array(
      const dynamic_array &other,
      parkway::utility::sort_order order =
        parkway::utility::sort_order::INCREASING) {
    this->sort_between_using_another_array(0, this->size() - 1, other, order);
  }

  inline void sort_between_using_another_array(
      std::size_t start, std::size_t end, const dynamic_array &other,
      parkway::utility::sort_order order =
        parkway::utility::sort_order::INCREASING) {
    parkway::utility::quick_sort_by_another_array(start, end, this->data(),
                                                  other.data(), order);
  }

  inline void random_permutation() {
    parkway::utility::random_permutation(this->data(), this->size());
  }

  inline int search(const_reference value) {
    for (std::size_t i = 0; i < this->size(); ++i) {
      if (data_->at(i) == value) {
        return i;
      }
    }
    return -1;
  }
};

}  // namespace data_structures
}  // namespace parkway

#endif  // DATA_STRUCTURES_DYNAMIC_ARRAY_HPP_
