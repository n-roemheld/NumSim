#include "Partitioning.h"

Partitioning::Partitioning()
{};

Partitioning::Partitioning(int MPI_rank, std::array<int,4> ranks_neighbors, std::array<bool,4> is_boundary, std::array<int,2> nCells, std::array<int,2> nCellsGlobal, std::array<int,2> nodeOffset) :
        MPI_rank_(MPI_rank), ranks_neighbors_(ranks_neighbors), is_boundary_(is_boundary), nCells_(nCells), nCellsGlobal_(nCellsGlobal), nodeOffset_(nodeOffset)
            {};

int Partitioning::ownRankNo()
{
    return MPI_rank_;
};

int Partitioning::rank_neighbor(int direction)
{
    return ranks_neighbors_[direction];
};

bool Partitioning::is_boundary(int direction) const
{
    return is_boundary_[direction];
};

std::array<int,2> Partitioning::nCells()
{
    return nCells_;
};
std::array<int,2> Partitioning::nCellsGlobal()
{
    return nCellsGlobal_;
};
bool Partitioning::ownPartitionContainsRightBoundary()
{
	return is_boundary(1);
}
bool Partitioning::ownPartitionContainsTopBoundary()
{
	return is_boundary(2);
}
std::array<int,2> Partitioning::nodeOffset()
{
	return nodeOffset_;
};
