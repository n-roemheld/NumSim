#include "GaussSeidelRedBlack.h"
#include <iostream>

GaussSeidelRedBlack::GaussSeidelRedBlack (std::shared_ptr< Discretization > discretization, double epsilon, int maximumNumberOfIterations) :
	PressureSolver(discretization, epsilon, maximumNumberOfIterations)
{};


void GaussSeidelRedBlack::solve()
{
	// std::cout << "GaussSeidelRedBlack" << std::endl;

	double dx = discretization_->dx();
	double dy = discretization_->dy();

	int it = 0;
	double res_squared = compute_res_parallel();

	while (it < maximumNumberOfIterations_ && res_squared > epsilon_*epsilon_)
	{
		// set boundary values for p to achieve 0-Neumann conditions
		setBoundaryValuesParallel();
		// perform one itertaion step on RED cells with (0,0) RED
		for(int j = discretization_->pJBegin(); j < discretization_->pJEnd(); j++)
		{
			for (int i = discretization_->pIBegin()+j%2; i < discretization_->pIEnd(); i+= 2)
			{
				discretization_->p(i,j) = (dx*dx*dy*dy)/(2*(dx*dx+dy*dy))
				* ( (discretization_->p(i-1,j) + discretization_->p(i+1,j)) / (dx*dx)
				+ (discretization_->p(i,j-1) + discretization_->p(i,j+1)) / (dy*dy)
				- discretization_->rhs(i,j) );
			};
		};
		pressure_communication();


    // perform one itertaion step on BLACK cells with (0,0) RED
    for(int j = discretization_->pJBegin(); j < discretization_->pJEnd(); j++)
    {
      for (int i = discretization_->pIBegin()+(j+1)%2; i < discretization_->pIEnd(); i+= 2)
      {
        discretization_->p(i,j) = (dx*dx*dy*dy)/(2*(dx*dx+dy*dy))
        * ( (discretization_->p(i-1,j) + discretization_->p(i+1,j)) / (dx*dx)
        + (discretization_->p(i,j-1) + discretization_->p(i,j+1)) / (dy*dy)
        - discretization_->rhs(i,j) );
      };
    };
		pressure_communication();


		// compute and update residuum
		res_squared = compute_res_parallel();
		it++;
	}
	setBoundaryValuesParallel();

};
