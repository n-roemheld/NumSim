#pragma once

#include "Discretization.h"

class CentralDifferences : public Discretization
{
public:

	CentralDifferences(std::array< int, 2 > nCells, std::array< double, 2 > meshWidth);

	virtual double computeDu2Dx (int i, int j) const;

	virtual double computeDv2Dy (int i, int j) const;

	virtual double computeDuvDx (int i, int j) const;

	virtual double computeDuvDy (int i, int j) const;

};
