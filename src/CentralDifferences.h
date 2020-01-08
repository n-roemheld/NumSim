#pragma once

#include "Discretization.h"

class CentralDifferences : public Discretization
{
public:

	// CentralDifferences(std::array< int, 2 > nCells, std::array< double, 2 > meshWidth);
    // CentralDifferences(std::array< int, 2 > nCells, std::array< double, 2 > meshWidth, std::shared_ptr<Array2D> geometryPVString, std::shared_ptr<Array2D> geometryPV1, std::shared_ptr<Array2D> geometryPV2, std::shared_ptr<Array2D> geometryTString, std::shared_ptr<Array2D> geometryT1);
    CentralDifferences(std::array< int, 2 > nCells, std::array< double, 2 > meshWidth, std::shared_ptr<Array2D> geometryPVString, std::shared_ptr<Array2D> geometryPVOrientation, std::shared_ptr<Array2D> geometryPV1, std::shared_ptr<Array2D> geometryPV2, std::shared_ptr<Array2D> geometryTString, std::shared_ptr<Array2D> geometryT1, Adapter& adapter, Settings& settings);

	virtual double computeDu2Dx (int i, int j) const;

	virtual double computeDv2Dy (int i, int j) const;

	virtual double computeDuvDx (int i, int j) const;

	virtual double computeDuvDy (int i, int j) const;

	virtual double computeDuTDx (int i, int j) const;

	virtual double computeDvTDy (int i, int j) const;

};
