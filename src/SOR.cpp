#include "SOR.h"
#include <iostream>

SOR::SOR (std::shared_ptr< Discretization > discretization, double epsilon, int maximumNumberOfIterations, double omega) :
PressureSolver(discretization, epsilon, maximumNumberOfIterations), omega_(omega)
{};

void SOR::solve()
{
	// change i,j for speed?
		std::array<double,2> mW = discretization_->meshWidth();
		double dx = mW[0];
		double dy = mW[1];

		int it = 0;
		// compute_res????
		// double res_squared = compute_res();
		double res_squared = 2*epsilon_*epsilon_;

		while (it <= maximumNumberOfIterations_ && res_squared > epsilon_*epsilon_)
		{
			// set boundary values for p to achieve 0-Neumann conditions
			setObstacleValues();
			setBoundaryValues();
			// perform one itertaion step
			for(int j = discretization_->pJBegin(); j < discretization_->pJEnd(); j++)
			{
				for (int i = discretization_->pIBegin(); i < discretization_->pIEnd(); i++)
				{
					// indices in geometry file (shifted by uIBegin and increased by one at the right (u) and upper(v) boundaries)
					int igeom = i-discretization_-> pIBegin()+1; // todo: double check!!
					int jgeom = j-discretization_-> pJBegin()+1;
					if(discretization_->geometryPVString(igeom, jgeom) == -1)
					{
						discretization_->p(i,j) = (1-omega_)*discretization_->p(i,j)+omega_*(dx*dx*dy*dy)/(2*(dx*dx+dy*dy))
						* ( (discretization_->p(i-1,j) + discretization_->p(i+1,j)) / (dx*dx)
						   + (discretization_->p(i,j-1) + discretization_->p(i,j+1)) / (dy*dy)
						   - discretization_->rhs(i,j) );
					}
				};
			};

			// compute and update residuum
			// wrong sign ????
			res_squared = compute_res();
			it++;
		}
		if(it > maximumNumberOfIterations_) std::cout << it << std::endl;
		setObstacleValues();
		setBoundaryValues();

};
