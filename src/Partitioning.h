#pragma once

#include <array>

class Partitioning
{
  public:
  Partitioning(int MPI_rank, std::array<int,4> ranks_neighbors, std::array<bool,4> is_boundary, std::array<int,2> nCells);
  int ownRankNo();
  bool is_boundary(int direction)// 0 = lower, 1 = right, 2 = upper, 3 = left
  std::array<int,4> ranks_neighbors(); // lower, right, upper, left
  // std::array<bool,4> is_boundary(); // lower, right, upper, left
  std::array<int,2> nCells();
  const std::array<int,2> nCellsGlobal();

  protected:
  const int MPI_rank_;
  const std::array<int,4> ranks_neighbors_; // lower, right, upper, left
  const std::array<bool,4> is_boundary_; // lower, right, upper, left
  const std::array<int,2> nCells_;
  const std::array<int,2> nCellsGlobal_;
};
