#include "ESGaussSeidel.h"

#include <iostream>

ESGaussSeidel::ESGaussSeidel (double epsilon, int maximumNumberOfIterations):
  EndSolver(epsilon, maximumNumberOfIterations)
{};

void ESGaussSeidel::solve(std::shared_ptr<MGGrid> mgg)
{
  std::array<double,2> mW = mgg->meshWidth();
  double dx = mW[0];
  double dy = mW[1];

  int it = 0;
	double res_squared = 2*epsilon_*epsilon_;

	while (it < maximumNumberOfIterations_ && res_squared > epsilon_*epsilon_)
	{
		// set boundary values for p to achieve 0-Neumann conditions
		setBoundaryValues(mgg);
		// perform one itertaion step
		for(int j = mgg->pJBegin(); j < mgg->pJEnd(); j++)
		{
			for (int i = mgg->pIBegin(); i < mgg->pIEnd(); i++)
			{
				mgg->p(i,j) = (dx*dx*dy*dy)/(2*(dx*dx+dy*dy))
				* ( (mgg->p(i-1,j) + mgg->p(i+1,j)) / (dx*dx)
				+ (mgg->p(i,j-1) + mgg->p(i,j+1)) / (dy*dy)
				- mgg->rhs(i,j) );
			};
		};
		// compute and update residuum
		res_squared = compute_res(mgg);
		it++;
	}
  // std::cout << "iteration " << it << std::endl;

	setBoundaryValues(mgg);

};
