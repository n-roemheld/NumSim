#include "GaussSeidel.h"

GaussSeidel::GaussSeidel (std::shared_ptr< Discretization > discretization, double epsilon, int maximumNumberOfIterations) :
	PressureSolver(discretization, epsilon, maximumNumberOfIterations)
{};


void GaussSeidel::solve()
{
	std::array<double,2> mW = discretization_->meshWidth();
	double dx = mW[0];
	double dy = mW[1];

	int it = 0;
	double res_squared = 2*epsilon_*epsilon_;

	while (it < maximumNumberOfIterations_ && res_squared > epsilon_*epsilon_)
	{
		// set boundary values for p to achieve 0-Neumann conditions
		setBoundaryValues();
		// perform one itertaion step
		for(int j = discretization_->pJBegin(); j < discretization_->pJEnd(); j++)
		{
			for (int i = discretization_->pIBegin(); i < discretization_->pIEnd(); i++)
			{
				discretization_->p(i,j) = (dx*dx*dy*dy)/(2*(dx*dx+dy*dy))
				* ( (discretization_->p(i-1,j) + discretization_->p(i+1,j)) / (dx*dx)
				+ (discretization_->p(i,j-1) + discretization_->p(i,j+1)) / (dy*dy)
				- discretization_->rhs(i,j) );
			};
		};
		// compute and update residuum
		res_squared = compute_res();
		it++;
	}
	setBoundaryValues();

};
