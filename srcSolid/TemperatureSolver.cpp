
#include "TemperatureSolver.h"
#include <math.h>
#include <memory>
#include <iostream>


TemperatureSolver::TemperatureSolver (std::shared_ptr< Discretization > discretization, double epsilon, int maximumNumberOfIterations) :
	discretization_(discretization), epsilon_(epsilon), maximumNumberOfIterations_(maximumNumberOfIterations)
{};

//!	set the boundary values to account for homogenous Neumann boundary conditions, this has to be called after every iteration
void TemperatureSolver::setBoundaryValues ()
{
	discretization_->setBoundaryValues_p(); 
};

// setObstacleValues is used, setObstacleValues2 not
// Marc
// void TemperatureSolver::setObstacleValues()
// {
// 	for(int j = discretization_->pJBegin(); j < discretization_->pJEnd(); j++)
// 	{
// 		for(int i = discretization_->pIBegin(); i < discretization_->pIEnd(); i++)
// 		{
// 			discretization_->setObstacleValues_p(i,j);
// 		}
// 	}
// }

// Nathanael
void TemperatureSolver::setObstacleValues2()
{
	for(int j = discretization_->pJBegin(); j < discretization_->pJEnd(); j++)
	{
		for(int i = discretization_->pIBegin(); i < discretization_->pIEnd(); i++)
		{
			// check if cell solid in setObstacleValues_p2
			discretization_->setObstacleValues_p2(i,j);
		}
	}
}

double TemperatureSolver::compute_res(double dt, double heatDiffusivity)
{
	//Array2D res_vec = Array2D(discretization_->nCells());
	double res = 0;
	FieldVariable T = discretization_->T();
	std::array<double,2> mW = discretization_->meshWidth();
	double dx = mW[0];
	double dy = mW[1];
	int cellCount = 0;
	//igeom and jgeom are identical to the indices of the pressure cell (?)
	for(int j = discretization_->pJBegin(); j < discretization_->pJEnd(); j++)
	{
		for (int i = discretization_->pIBegin(); i < discretization_->pIEnd(); i++)
		{
			if(discretization_->geometryPVString(i,j) == -1) // proof, if cell is fluid cell
			{
				res += pow(T(i,j) - dt * heatDiffusivity * ((T(i+1,j)-2*T(i,j)+T(i-1,j))/(dx*dx) + (T(i,j+1)-2*T(i,j)+T(i,j-1))/(dy*dy)) - discretization_->rhs(i,j),2);
				cellCount++;
			}
		};
	};
	// average residuum with respect to number of cells
	res /= cellCount;

	return res;
};
