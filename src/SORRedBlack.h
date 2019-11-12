#pragma once

#include "PressureSolver.h"

class SORRedBlack : public PressureSolver
{
public:
	SORRedBlack (std::shared_ptr< Discretization > discretization, double epsilon, int maximumNumberOfIterations, double omega);

	//! solve the system of the Poisson equation for pressure
	void solve();

private:
	double omega_;
};
