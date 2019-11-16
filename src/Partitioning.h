#pragma once

#include <array>

class Partitioning
{
  public:
  Partitioning();
  Partitioning(int MPI_rank, std::array<int,4> ranks_neighbors, std::array<bool,4> is_boundary, std::array<int,2> nCells, std::array<int,2> nCellsGlobal);
  int ownRankNo();
  bool is_boundary(int direction) const;// 0 = lower, 1 = right, 2 = upper, 3 = left
  int rank_neighbor(int direction);
  // std::array<int,4> ranks_neighbors(); // lower, right, upper, left
  // std::array<bool,4> is_boundary(); // lower, right, upper, left
  std::array<int,2> nCells();
  //const
  std::array<int,2> nCellsGlobal();

  protected:
    // eig const
  int MPI_rank_;
  std::array<int,4> ranks_neighbors_; // lower, right, upper, left
  std::array<bool,4> is_boundary_; // lower, right, upper, left
  std::array<int,2> nCells_;
  std::array<int,2> nCellsGlobal_;
};
