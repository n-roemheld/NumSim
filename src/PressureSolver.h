#pragma once

#include "Discretization.h"
#include <memory>

class PressureSolver
{
public:
	PressureSolver (std::shared_ptr<Discretization> discretization, double epsilon, int maximumNumberOfIterations);

	//!solve the Poisson problem for the pressure, using the rhs and p field variables in the staggeredGrid
	virtual void solve () = 0;

protected:
	//!	set the boundary values to account for homogenous Neumann boundary conditions, this has to be called after every iteration
	void setBoundaryValues ();

	std::shared_ptr< Discretization > discretization_;
	double 	epsilon_;
	int 	maximumNumberOfIterations_;

	double compute_res();
};
