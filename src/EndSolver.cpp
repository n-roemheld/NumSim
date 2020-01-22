# pragma once 

#include "EndSolver.h"
#include <math.h>

EndSolver::EndSolver(double epsilon, int maximumNumberOfIterations) : 
    epsilon_(epsilon), maximumNumberOfIterations_(maximumNumberOfIterations)#
{};

void EndSolver::setBoundaryValues(std::shared_ptr<MGGrid> mgg)
{
    // lower p ghost layer
		int j = mgg->pJBegin()-1;
		for(int i = mgg->pIBegin(); i < mgg->pIEnd(); i++)
		{
			mgg->p(i,j) = mgg->p(i,j+1);
		};
		// upper p ghost layer
		j = mgg->pJEnd();
		for(int i = mgg->pIBegin(); i < mgg->pIEnd(); i++)
		{
			mgg->p(i,j) = mgg->p(i,j-1);
		};
		// left p ghost layer
		int i = mgg->pIBegin()-1;
		for(int j = mgg->pJBegin()-1; j < mgg->pJEnd()+1; j++)
		{
			mgg->p(i,j) = mgg->p(i+1,j);
		}
		// right p ghost layer
		i = mgg->pIEnd();
		for(int j = mgg->pJBegin()-1; j < mgg->pJEnd()+1; j++)
		{
			mgg->p(i,j) = mgg->p(i-1,j);
		}
};

double EndSolver::compute_res(std::shared_ptr<MGGrid> mgg)
{
    //Array2D res_vec = Array2D(mgg->nCells());
	double res = 0;
	FieldVariable p = mgg->p();
	std::array<double,2> mW = mgg->meshWidth();
	double dx = mW[0];
	double dy = mW[1];
	for(int j = mgg->pJBegin(); j < mgg->pJEnd(); j++)
	{
		for (int i = mgg->pIBegin(); i < mgg->pIEnd(); i++)
		{
			res += pow((p(i+1,j)-2*p(i,j)+p(i-1,j))/(dx*dx) + (p(i,j+1)-2*p(i,j)+p(i,j-1))/(dy*dy) - mgg->rhs(i,j),2);
		};
	};

	// average residuum with respect to number of cells
	res /= (mgg->nCells()[0]*mgg->nCells()[1]);
	return res;
};