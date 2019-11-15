#include "Partitioning.h"

Partitioning::Partitioning(int MPI_rank, std::array<int,4> ranks_neighbors, std::array<bool,4> is_boundary, std::array<int,2> nCells) :
        MPI_rank_(MPI_rank), ranks_neighbors_(ranks_neighbors), is_boundary_(is_boundary), nCells_(nCells) {};

int Partitioning::ownRankNo()
{
    return MPI_rank_;
};

int Partitioning::rank_neighbors(int direction)
{
    return ranks_neighbors_[direction];
};

bool Partitioning::is_boundary(int direction)
{
    return is_boundary_[direction];
};
std::array<int,2> Partitioning::nCells()
{
    return nCells_;
};
