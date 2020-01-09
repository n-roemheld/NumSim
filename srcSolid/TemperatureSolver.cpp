
#include "TemperatureSolver.h"
#include <math.h>
#include <memory>
#include <iostream>


TemperatureSolver::TemperatureSolver (std::shared_ptr< Discretization > discretization, double epsilon, int maximumNumberOfIterations) :
	discretization_(discretization), epsilon_(epsilon), maximumNumberOfIterations_(maximumNumberOfIterations)
{};

//!	set the boundary values to account for homogenous Neumann boundary conditions, this has to be called after every iteration
void TemperatureSolver::setBoundaryValues (std::vector<double> readData)
{
	// discretization_->setBoundaryValues_p(); 
	// todo: double check indices

	// locations:
	int left = 0;
	int right = 1;
	int lower = 2;
	int upper = 3;

	// setting T boundaries without corners (p grid, all ghost cells)
	// lower T
	int j = discretization_->pJBegin()-1;
	for(int i = discretization_->pIBegin(); i < discretization_->pIEnd(); i++)
	{
		discretization_->setBoundaryValues_T(lower,i,j);
	};
	// upper T
	j = discretization_->pJEnd();
	for(int i = discretization_->pIBegin(); i < discretization_->pIEnd(); i++)
	{
		discretization_->setBoundaryValues_T(upper,i,j);
	};
	// left T
	int i = discretization_->pIBegin()-1;
	for(int j = discretization_->pJBegin()-1; j < discretization_->pJEnd()+1; j++)
	{
		discretization_->setBoundaryValues_T(left,i,j);
	}
	// right T
	i = discretization_->pIEnd();
	for(int j = discretization_->pJBegin()-1; j < discretization_->pJEnd()+1; j++)
	{
		discretization_->setBoundaryValues_T(right,i,j);
	}

	for(int v = 0; v < discretization_->settings_->vertexSize; v++)
	{
		int i = discretization_->settings_->vertex_i.at(v);
		int j = discretization_->settings_->vertex_j.at(v);

		// std::cout << "here!" << std::endl;

		std::array<double,2> mW = discretization_->meshWidth();
		double dx = mW[0];
		double dy = mW[1];

		// neighbour indices
		int in = i;
		int jn = j;
		double h = 0;

		// std::cout << "c" << std::endl;
		// std::cout << "osize" << discretization_->settings_->orientation_.size() << std::endl;


		int orientation = discretization_->settings_->orientation_.at(v);

		// std::cout << "here2" << std::endl;

		switch (orientation)
		{
			case 0: std::cout << "no orientation assigned!" << std::endl; break;
			case 1: in = i-1; h = dx; break;
			case 2: in = i+1; h = dx; break;
			case 3: jn = j-1; h = dy; break;
			case 4: jn = j+1; h = dy; break;
			case 5:
			case 6:
			case 7:
			case 8: h = 0; break;
			default: std::cout << "unknown orientation" << std::endl; break;
		}

		// std::cout << "here3!" << std::endl;


		discretization_->T(i,j) = h*readData.at(v)*discretization_->settings_->re*discretization_->settings_->prandtl        //lkjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjj
		                        + discretization_->T(in,jn);
	}
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
