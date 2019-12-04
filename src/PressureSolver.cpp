
#include "PressureSolver.h"
#include <math.h>
#include <memory>


PressureSolver::PressureSolver (std::shared_ptr< Discretization > discretization, double epsilon, int maximumNumberOfIterations) :
	discretization_(discretization), epsilon_(epsilon), maximumNumberOfIterations_(maximumNumberOfIterations)
{};

//!	set the boundary values to account for homogenous Neumann boundary conditions, this has to be called after every iteration
void PressureSolver::setBoundaryValues ()
{
	// locations:
	int left = 0;
	int right = 1;
	int lower = 2;
	int upper = 3;

	// setting p boundaries without corners (p grid, all ghost cells)
	// lower p
	int j = discretization_->pJBegin()-1;
	for(int i = discretization_->pIBegin(); i < discretization_->pIEnd(); i++)
	{
		discretization_->setBoundaryValues_p(lower,i,j);
	};
	// upper p
	j = discretization_->pJEnd();
	for(int i = discretization_->pIBegin(); i < discretization_->pIEnd(); i++)
	{
		discretization_->setBoundaryValues_p(upper,i,j);
	};
	// left p
	int i = discretization_->pIBegin()-1;
	for(int j = discretization_->pJBegin()-1; j < discretization_->pJEnd()+1; j++)
	{
		discretization_->setBoundaryValues_p(left,i,j);
	}
	// right p
	i = discretization_->pIEnd();
	for(int j = discretization_->pJBegin()-1; j < discretization_->pJEnd()+1; j++)
	{
		discretization_->setBoundaryValues_p(right,i,j);
	}
};

void PressureSolver::setObstacleValues()
{
	for(int j = discretization_->pJBegin(); j < discretization_->pJEnd(); j++)
	{
		for(int i = discretization_->pIBegin(); i < discretization_->pIEnd(); i++)
		{
			discretization_->setObstacleValues_p(i,j);
		}
	}
}

void PressureSolver::setObstacleValues2()
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

double PressureSolver::compute_res()
{
	//Array2D res_vec = Array2D(discretization_->nCells());
	double res = 0;
	FieldVariable p = discretization_->p();
	std::array<double,2> mW = discretization_->meshWidth();
	double dx = mW[0];
	double dy = mW[1];
	int fluidCellCount = 0;
	//igeom and jgeom are identical to the indices of the pressure cell (?)
	for(int j = discretization_->pJBegin(); j < discretization_->pJEnd(); j++)
	{
		for (int i = discretization_->pIBegin(); i < discretization_->pIEnd(); i++)
		{
			if(discretization_->geometryPVString(i,j) == -1) // proof, if cell is fluid cell
			{
				res += pow((p(i+1,j)-2*p(i,j)+p(i-1,j))/(dx*dx) + (p(i,j+1)-2*p(i,j)+p(i,j-1))/(dy*dy) - discretization_->rhs(i,j),2);
				fluidCellCount++;
			}
		};
	};

	// average residuum with respect to number of cells
	res /= fluidCellCount;
	return res;
};
