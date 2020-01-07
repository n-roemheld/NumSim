#pragma once

#include "TemperatureSolver.h"

class SOR : public TemperatureSolver
{
public:
	SOR (std::shared_ptr< Discretization > discretization, double epsilon, int maximumNumberOfIterations, double omega);

	//! solve the system of the Poisson equation for pressure
	void solve();

private:
	double omega_;
};
