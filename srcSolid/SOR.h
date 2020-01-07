#pragma once

#include "PressureSolver.h"

class SOR : public PressureSolver
{
public:
	SOR (std::shared_ptr< Discretization > discretization, double epsilon, int maximumNumberOfIterations, double omega);

	//! solve the system of the Poisson equation for pressure
	void solve();

private:
	double omega_;
};
