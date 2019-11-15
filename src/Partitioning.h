#pragma once

#include <array>

class Partitioning
{
  public:
  Partitioning(int MPI_rank, std::array<int,4> ranks_neighbors, std::array<bool,4> is_boundary, std::array<int,2> nCells);
  int rank();
  std::array<int,4> ranks_neighbors();
  std::array<bool,4> is_boundary();
  std::array<int,2> nCells();

  protected:
  const int MPI_rank_;
  const std::array<int,4> ranks_neighbors_; // lower, right, upper, left
  const std::array<bool,4> is_boundary_; // lower, right, upper, left
  const std::array<int,2> nCells_;
};
