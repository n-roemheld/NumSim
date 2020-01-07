#include "GaussSeidel.h"

GaussSeidel::GaussSeidel (std::shared_ptr< Discretization > discretization, double epsilon, int maximumNumberOfIterations) :
	TemperatureSolver(discretization, epsilon, maximumNumberOfIterations)
{};


void GaussSeidel::solve(double dt, double heatDiffusivity)
{
	std::array<double,2> mW = discretization_->meshWidth();
	double dx = mW[0];
	double dy = mW[1];

	int it = 0;
	setObstacleValues2();
	setBoundaryValues();
	double res_squared = compute_res(dt, heatDiffusivity);
	// double res_squared = 2*epsilon_*epsilon_;

	while (it < maximumNumberOfIterations_ && res_squared > epsilon_*epsilon_)
	{
		// set boundary values for p to achieve 0-Neumann conditions
		setObstacleValues2();
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
					discretization_->T(i,j) = 1./(1 + dt * heatDiffusivity * (2./(dx * dx) + 2./(dy * dy))) 
												* (discretization_->rhs(i,j) + dt * heatDiffusivity 
												* ((discretization_->T(i+1,j) + discretization_->T(i-1,j))/(dx*dx) 
												+ (discretization_->T(i,j+1) + discretization_->T(i,j-1))/(dy*dy)));
				}
			};
		};
		// compute and update residuum
		res_squared = compute_res(dt, heatDiffusivity);
		it++;
	}
	setObstacleValues2();
	setBoundaryValues();

};
