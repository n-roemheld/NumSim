#include "Array2D.h"

#include <cassert>
#include <assert.h>
#include <iostream>

Array2D::Array2D(std::array<int,2> size) :
  size_(size)
{
  // allocate data, initialize to 0
  data_.resize(size_[0]*size_[1], 0.0);
};

//! get the size
std::array<int,2> Array2D::size() const
{
  return size_;
};

double &Array2D::operator()(int i, int j)
{
  const int index = j*size_[0] + i;
  // std::cout << "i" << i  << "j" << j << std::endl;
  // std::cout << size_[0] << size_[1] << std::endl;
  // assert that indices are in range
  if (i<0 || j<0)
  {
    std::cout << "i: " << i << ", j: " << j << std::endl;
  }
  if (i>=size_[0] || j>=size_[1] )
  {
    std::cout << "i: " << i << ", j: " << ", size0: " << size_[0] << ", size1: " << size_[1] << std::endl;
  }
  assert(0 <= i && i < size_[0]);
  assert(0 <= j && j < size_[1]);
  assert(j*size_[0] + i < (int)data_.size());
  assert(0 <= i && i < size_[0]);
  assert(0 <= j && j < size_[1]);
  assert(j*size_[0] + i < (int)data_.size());

  return data_[index];
};

double Array2D::operator()(int i, int j) const
{
  const int index = j*size_[0] + i;
  // std::cout << "i" << i  << "j" << j << std::endl;
  // std::cout << size_[0] << size_[1] << std::endl;
  // assert that indices are in range
  if (i<0 || j<0)
  {
    std::cout << "const: i: " << i << ", j: " << j << std::endl;
  }
  if (i>=size_[0] || j>=size_[1] )
  {
    std::cout << "const: i: " << i << ", j: " << ", size0: " << size_[0] << ", size1: " << size_[1] << std::endl;
  }
  assert(0 <= i && i < size_[0]);
  assert(0 <= j && j < size_[1]);
  assert(j*size_[0] + i < (int)data_.size());

  return data_[index];
};

double *Array2D::data() 
{
  return data_.data();
}
