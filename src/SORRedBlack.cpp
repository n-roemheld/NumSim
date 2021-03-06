#include "SORRedBlack.h"
#include <iostream>

SORRedBlack::SORRedBlack (std::shared_ptr< Discretization > discretization, double epsilon, int maximumNumberOfIterations, double omega) :
PressureSolver(discretization, epsilon, maximumNumberOfIterations), omega_(omega)
{};

void SORRedBlack::solve()
{
	// std::cout << "SORRedBlack" << std::endl;
	// change i,j for speed?
	double dx = discretization_->dx();
	double dy = discretization_->dy();

		int it = 0;
		// compute_res????
		double res_squared = compute_res_parallel();
		// double res_squared = 2*epsilon_*epsilon_;

		while (it <= maximumNumberOfIterations_ && res_squared > epsilon_*epsilon_)
		{
			// set boundary values for p to achieve 0-Neumann conditions
			setBoundaryValuesParallel();
			// perform one itertaion step for RED with (0,0) RED
			for(int j = discretization_->pJBegin(); j < discretization_->pJEnd(); j++)
			{
				for (int i = discretization_->pIBegin()+j%2; i < discretization_->pIEnd(); i+= 2)
				{
					discretization_->p(i,j) = (1-omega_)*discretization_->p(i,j)+omega_*(dx*dx*dy*dy)/(2*(dx*dx+dy*dy))
					* ( (discretization_->p(i-1,j) + discretization_->p(i+1,j)) / (dx*dx)
					   + (discretization_->p(i,j-1) + discretization_->p(i,j+1)) / (dy*dy)
					   - discretization_->rhs(i,j) );
				};
			};
            pressure_communication();

			// perform one itertaion step for BLACK with (0,0) RED
			for(int j = discretization_->pJBegin(); j < discretization_->pJEnd(); j++)
			{
				for (int i = discretization_->pIBegin()+(j+1)%2; i < discretization_->pIEnd(); i+= 2)
				{
					discretization_->p(i,j) = (1-omega_)*discretization_->p(i,j)+omega_*(dx*dx*dy*dy)/(2*(dx*dx+dy*dy))
					* ( (discretization_->p(i-1,j) + discretization_->p(i+1,j)) / (dx*dx)
						 + (discretization_->p(i,j-1) + discretization_->p(i,j+1)) / (dy*dy)
						 - discretization_->rhs(i,j) );
				};
			};
            pressure_communication();

			// compute and update residuum
			// wrong sign ????
			res_squared = compute_res_parallel();
			it++;
		}
		if(it > maximumNumberOfIterations_) std::cout << it << std::endl;
		setBoundaryValuesParallel();

};
