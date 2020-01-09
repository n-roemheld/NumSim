#pragma once

#include <array>
#include "StaggeredGrid.h"


class Discretization : public StaggeredGrid
{
public:
	// Discretization (std::array< int, 2 > nCells, std::array< double, 2 > meshWidth);
    // Discretization (std::array< int, 2 > nCells, std::array< double, 2 > meshWidth, std::shared_ptr<Array2D> geometryPVString, std::shared_ptr<Array2D> geometryPV1, std::shared_ptr<Array2D> geometryPV2, std::shared_ptr<Array2D> geometryTString, std::shared_ptr<Array2D> geometryT1);
    Discretization (std::array< int, 2 > nCells, std::array< double, 2 > meshWidth, std::shared_ptr<Array2D> geometryPVString, std::shared_ptr<Array2D> geometryPVOrientation, std::shared_ptr<Array2D> geometryPV1, std::shared_ptr<Array2D> geometryPV2, std::shared_ptr<Array2D> geometryTString, std::shared_ptr<Array2D> geometryT1, std::shared_ptr<Settings> settings);

	virtual double computeDu2Dx (int i, int j) const = 0;

	virtual double computeDv2Dy (int i, int j) const = 0;

	virtual double computeDuvDx (int i, int j) const = 0;

	virtual double computeDuvDy (int i, int j) const = 0;

	virtual double computeD2uDx2 (int i, int j) const;

	virtual double computeD2uDy2 (int i, int j) const;

	virtual double computeD2vDx2 (int i, int j) const;

	virtual double computeD2vDy2 (int i, int j) const;

	virtual double computeDpDx (int i, int j) const;

	virtual double computeDpDy (int i, int j) const;

	double computeD2TDx2 (int i, int j) const;

	double computeD2TDy2 (int i, int j) const;

	virtual double computeDuTDx (int i, int j) const=0;

	virtual double computeDvTDy (int i, int j) const=0;
};
