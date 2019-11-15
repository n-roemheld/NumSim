#include "Partitioning.h"

Partitioning::Partitioning(int MPI_rank, std::array<int,4> ranks_neighbors, std::array<bool,4> is_boundary, std::array<int,2> nCells) :
        MPI_rank_(MPI_rank), ranks_neighbors_(ranks_neighbors), is_boundary_(is_boundary), nCells_(nCells) {};

int rank() {
    return rank_;
};
std::array<int,4> Partitioning::ranks_neighbors()
{
    return ranks_neighbors;
};
std::array<bool,4> Partitioning::is_boundary()
{
    return is_boundary_;
};
std::array<int,2> Partitioning::nCells()
{
    return nCells_;
};
