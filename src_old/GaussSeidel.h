#pragma once

#include "PressureSolver.h"

class GaussSeidel : public PressureSolver
{
public:

	GaussSeidel (std::shared_ptr< Discretization > discretization, double epsilon, int maximumNumberOfIterations);


	void solve ();

};
