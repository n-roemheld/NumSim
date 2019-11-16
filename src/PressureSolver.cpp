
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

void PressureSolver::setBoundaryValuesParallel()
{
    // neigbor indices
    int below = 0;
    int above = 2;
    int right = 1;
    int left = 3;

    // lower p ghost layer
    if (discretization_->is_boundary(below))
    {
        int j = discretization_->pJBegin()-1;
        for(int i = discretization_->pIBegin(); i < discretization_->pIEnd(); i++)
        {
            discretization_->p(i,j) = discretization_->p(i,j+1);
        };
    }

// upper p ghost layer
    if (discretization_->is_boundary(above))
    {
        int j = discretization_->pJEnd();
        for(int i = discretization_->pIBegin(); i < discretization_->pIEnd(); i++)
        {
            discretization_->p(i,j) = discretization_->p(i,j-1);
        };
    }

// left p ghost layer
    if (discretization_->is_boundary(left))
    {
        int i = discretization_->pIBegin()-1;
        for(int j = discretization_->pJBegin()-1; j < discretization_->pJEnd()+1; j++)
        {
            discretization_->p(i,j) = discretization_->p(i+1,j);
        }
    }

// right p ghost layer
    if (discretization_->is_boundary(right))
    {
        int i = discretization_->pIEnd();
        for(int j = discretization_->pJBegin()-1; j < discretization_->pJEnd()+1; j++)
        {
            discretization_->p(i,j) = discretization_->p(i-1,j);
        }
    }

};

void PressureSolver::pressure_communication()
{
    // neigbor indices
    int below = 0;
    int above = 2;
    int right = 1;
    int left = 3;

    // MPI tags of information sent (communication directions)
    int down = 0;
    int to_right = 1;
    int up = 2;
    int to_left = 3;

    // Sending Pressures
    // send to right
    discretization_->send_boundary_horizontal_p(to_right, discretization_->pIEnd()-1, discretization_->pJBegin(), discretization_->pJEnd(), discretization_->rank_neighbor(right), discretization_->is_boundary(right));
    // send to left
    discretization_->send_boundary_horizontal_p(to_left, discretization_->pIBegin(), discretization_->pJBegin(), discretization_->pJEnd(), discretization_->rank_neighbor(left), discretization_->is_boundary(left));
    // send to above
    discretization_->send_boundary_vertical_p(up, discretization_->pJEnd()-1, discretization_->pIBegin(), discretization_->pIEnd(), discretization_->rank_neighbor(above), discretization_->is_boundary(above));
    // send to below
    discretization_->send_boundary_vertical_p(down, discretization_->pJBegin(), discretization_->pIBegin(), discretization_->pIEnd(), discretization_->rank_neighbor(below), discretization_->is_boundary(below));

    // Receiving Pressures
    std::vector<MPI_Request> requests;
		MPI_Request current_request;
    // receive from right
    current_request = discretization_->receive_boundary_horizontal_p(to_left, discretization_->pIEnd(), discretization_->pJBegin(), discretization_->pJEnd(), discretization_->rank_neighbor(right), discretization_->is_boundary(right));
    requests.push_back(current_request);
    // receive from left
    current_request = discretization_->receive_boundary_horizontal_p(to_right, discretization_->pIBegin()-1, discretization_->pJBegin(), discretization_->pJEnd(), discretization_->rank_neighbor(left), discretization_->is_boundary(left));
    requests.push_back(current_request);
    // receive from above
    current_request = discretization_->receive_boundary_vertical_p(down, discretization_->pJEnd(), discretization_->pIBegin(), discretization_->pIEnd(), discretization_->rank_neighbor(above), discretization_->is_boundary(above));
    requests.push_back(current_request);
    // receive from below
    current_request = discretization_->receive_boundary_vertical_p(up, discretization_->pJBegin()-1, discretization_->pIBegin(), discretization_->pIEnd(), discretization_->rank_neighbor(below), discretization_->is_boundary(below));
    requests.push_back(current_request);

		// std::vector<MPI_Status> statuse(4,MPI_STATUS_IGNORE);
    // MPI_Waitall(4, &requests, statuse);
		MPI_Waitall(requests.size(), requests.data(), MPI_STATUS_IGNORE);


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

double PressureSolver::compute_res_parallel()
{
    double res_local = compute_res()*(discretization_->nCells()[0]*discretization_->nCells()[1]);
    double res_global;

    MPI_Allreduce(&res_local, &res_global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    res_global /= (discretization_->nCellsGlobal()[0]*discretization_->nCellsGlobal()[1]);
    return res_global;
}
