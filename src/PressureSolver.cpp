
#include "PressureSolver.h"
#include <math.h>
#include <memory>


PressureSolver::PressureSolver (std::shared_ptr< Discretization > discretization, double epsilon, int maximumNumberOfIterations) :
	discretization_(discretization), epsilon_(epsilon), maximumNumberOfIterations_(maximumNumberOfIterations)
{};

//!	set the boundary values to account for homogenous Neumann boundary conditions, this has to be called after every iteration
void PressureSolver::setBoundaryValues ()
{
	// ???????? verschattet C++ Variablen??? -> scheinbar ja
	// p setting
		// lower p ghost layer
		int j = discretization_->pJBegin()-1;
		for(int i = discretization_->pIBegin(); i < discretization_->pIEnd(); i++)
		{
			discretization_->p(i,j) = discretization_->p(i,j+1);
		};
		// upper p ghost layer
		j = discretization_->pJEnd();
		for(int i = discretization_->pIBegin(); i < discretization_->pIEnd(); i++)
		{
			discretization_->p(i,j) = discretization_->p(i,j-1);
		};
		// left p ghost layer
		int i = discretization_->pIBegin()-1;
		for(int j = discretization_->pJBegin()-1; j < discretization_->pJEnd()+1; j++)
		{
			discretization_->p(i,j) = discretization_->p(i+1,j);
		}
		// right p ghost layer
		i = discretization_->pIEnd();
		for(int j = discretization_->pJBegin()-1; j < discretization_->pJEnd()+1; j++)
		{
			discretization_->p(i,j) = discretization_->p(i-1,j);
		}
};

double PressureSolver::compute_res()
{
	//Array2D res_vec = Array2D(discretization_->nCells());
	double res = 0;
	FieldVariable p = discretization_->p();
	std::array<double,2> mW = discretization_->meshWidth();
	double dx = mW[0];
	double dy = mW[1];
	for(int j = discretization_->pJBegin(); j < discretization_->pJEnd(); j++)
	{
		for (int i = discretization_->pIBegin(); i < discretization_->pIEnd(); i++)
		{
			res += pow((p(i+1,j)-2*p(i,j)+p(i-1,j))/(dx*dx) + (p(i,j+1)-2*p(i,j)+p(i,j-1))/(dy*dy) - discretization_->rhs(i,j),2);
		};
	};

	// average residuum with respect to number of cells
	res /= (discretization_->nCells()[0]*discretization_->nCells()[1]);
	return res;
};
