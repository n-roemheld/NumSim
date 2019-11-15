#pragma once

#include "PressureSolver.h"

class GaussSeidelRedBlack : public PressureSolver
{
public:

	GaussSeidelRedBlack (std::shared_ptr< Discretization > discretization, double epsilon, int maximumNumberOfIterations);


	void solve ();

};
