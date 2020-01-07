#pragma once

#include "TemperatureSolver.h"

class GaussSeidel : public TemperatureSolver
{
public:

	GaussSeidel (std::shared_ptr< Discretization > discretization, double epsilon, int maximumNumberOfIterations);


	void solve (double dt, double heatDiffusivity);

};
