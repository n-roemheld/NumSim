#pragma once

#include "DonorCell.h"

class DonorCellParallel : public DonorCell
{
public:
	DonorCellParallel (std::array< int, 2 > nCells, std::array< double, 2 > meshWidth, double alpha, Partitioning parti);

private:
	double alpha_;
};
