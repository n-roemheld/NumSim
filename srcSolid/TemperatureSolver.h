#pragma once

#include "Discretization.h"
#include <memory>

class TemperatureSolver
{
public:
	TemperatureSolver (std::shared_ptr<Discretization> discretization, double epsilon, int maximumNumberOfIterations);

	//!solve the Poisson problem for the pressure, using the rhs and p field variables in the staggeredGrid
	virtual void solve (double dt, double heatDiffusivity) = 0;

protected:
	//!	set the boundary values to account for homogenous Neumann boundary conditions, this has to be called after every iteration
	void setBoundaryValues ();

	// void setObstacleValues();
	void setObstacleValues2();

	std::shared_ptr< Discretization > discretization_;
	double 	epsilon_;
	int 	maximumNumberOfIterations_;

	double compute_res(double dt, double heatDiffusivity);
};
