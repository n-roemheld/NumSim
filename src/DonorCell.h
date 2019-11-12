#pragma once

#include "Discretization.h"

class DonorCell : public Discretization
{
public:
	DonorCell (std::array< int, 2 > nCells, std::array< double, 2 > meshWidth, double alpha);

	virtual double computeDu2Dx (int i, int j) const;

	virtual double computeDv2Dy (int i, int j) const;

	virtual double computeDuvDx (int i, int j) const;

	virtual double computeDuvDy (int i, int j) const;

private:
	double alpha_;
};
