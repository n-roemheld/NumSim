
#include "ComputationParallel.h"
#include "Partitioning.h"

#include <math.h>
#include <cassert>
#include <assert.h>
#include <vector>

void ComputationParallel::initialize (int argc, char *argv[])
{
    //Computation::initialize (int argc, char *argv[]);
    settings_.loadFromFile(argv[1]);
	// computing meshWidth (everywhere with global values)
	double dx = settings_.physicalSize[0]/settings_.nCells[0];
	double dy = settings_.physicalSize[1]/settings_.nCells[1];
	meshWidth_ = {dx, dy};

    
    MPI_Init(Null,Null); // Number of processes determined by command line
    // Get the number of processes
    int MPI_n_processes;
    MPI_Comm_size(MPI_COMM_WORLD, &MPI_n_processes);
    // Get the rank of the process
    int MPI_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &MPI_rank);
    
    // compute nCells and assign physical relationships
    if(MPI_rank == 0)
    {
        settings_.printSettings();
        // Number of partitions in both dimensions
        int n_pars_x = int (std::sqrt(settings_->nCells[0]));
        int n_pars_y = int (std::sqrt(settings_->nCells[1]));
        assert( n_pars_x*n_pars_y == MPI_n_processes);
        // Number of cells in each partiotion except the last one
        int n_Cells_sub_x = int (settings_->nCells[0]) / n_pars_x);
        int n_Cells_sub_y = int (settings_->nCells[1]) / n_pars_y);
        // Number of cells in the last partitions in each dimension 
        int n_Cells_sub_x_last = settings_->nCells[0] - n_pars_x*n_Cells_sub_x;
        int n_Cells_sub_y_last = settings_->nCells[1] - n_pars_y*n_Cells_sub_y;
        
        // Defining ranks for all partitions
        Array2D ranks_domain({n_pars_x,n_pars_y});
        int rank_counter = 0;
        for(int j = 0; j < n_pars_y; j++)
        {
            for (int i = 0; i < n_pars_x; i++)
            {
                ranks_domain(i,j) = rank_counter;
                rank_counter++;
            }
        }
        
        // Setting neighbor and boundary properties for all partitions
        for(int j = 0; j < n_pars_y; j++)
        {
            for (int i = 0; i < n_pars_x; i++)
            {
                if (i == 0 && j == 0)
                {
                    std::vector<int> ranks_neighbors{ranks_domain(i,j-1), ranks_domain(i+1,j), ranks_domain(i,j+1), ranks_domain(i-1,j)}; // bottom, right, upper, left; caution: check for limits (boundaries)!
                    std::vector<int> is_boundary{ (j-1) == -1, (i+1) == n_pars_x, (j+1) == n_pars_y, (i-1) == -1};
                    std::vector<int> nCells_sub{n_Cells_sub_x,n_Cells_sub_y};
                    // Not changing physicalSize because it's not used. Caution: Inconsistent to nCells
                    settings_.nCells = nCells_sub;
                    discretization_->set_partitioning(MPI_rank, ranks_neighbors, is_boundary, nCells_sub, settings_->nCells);
                }
                else
                {
                    std::vector<int> ranks_neighbors{ranks_domain(i,j-1), ranks_domain(i+1,j), ranks_domain(i,j+1), ranks_domain(i-1,j)}; // bottom, right, upper, left; caution: check for limits (boundaries)!
                    std::vector<int> is_boundary{ (j-1) == -1, (i+1) == n_pars_x, (j+1) == n_pars_y, (i-1) == -1};
                    std::vector<int> nCells_sub{0,0};
                    if (j == n_pars_y) {
                        nCells_sub[1] = n_Cells_sub_y_last;
                    } else {
                        nCells_sub[1] = n_Cells_sub_y;
                    }
                    if (i == n_pars_x) {
                        nCells_sub[0] = n_Cells_sub_x_last;
                    } else {
                        nCells_sub[0] = n_Cells_sub_x;
                    }
                    
                    // send part to partition ranks_domain(i,j) or store in arrays and broadcast
                    MPI_Isend(&ranks_neighbors, 10, MPI_INT, ranks_domain(i,j),0,MPI_COMM_WORLD);
                    MPI_Isend(&is_boundary, 10, MPI_INT, ranks_domain(i,j),1,MPI_COMM_WORLD);
                    MPI_Isend(&nCells_sub, 10, MPI_INT, ranks_domain(i,j),2,MPI_COMM_WORLD);
                }
            }
        }
    
    }
    else
    {
        std::vector<int> ranks_neighbors; // bottom, right, upper, left; caution: check for limits (boundaries)!
        std::vector<int> is_boundary;
        std::vector<int> nCells_sub;
        std::vector<MPI_Request> requests;
        MPI_Request current_request;
        MPI_Irecv(&ranks_neighbors, 10, MPI_INT, 0, 0, MPI_COMM_WORLD, &current_request);
        requests.push_back(current_request);
        MPI_Irecv(&is_boundary, 10, MPI_INT, 0, 1, MPI_COMM_WORLD, &current_request);
        requests.push_back(current_request);
        MPI_Irecv(&nCells_sub, 10, MPI_INT, 0, 2, MPI_COMM_WORLD, &current_request);
        requests.push_back(current_request);
        discretization_->set_partitioning(MPI_rank, ranks_neighbors, is_boundary, nCells_sub, settings_->nCells);
        MPI_Wait(&requests,MPI_STATUS_IGNORE);
        // Not changing physicalSize because it's not used. Caution: Inconsistent to nCells
        settings_.nCells = nCells_sub;
    }
    
};

void ComputationParallel::computeTimeStepWidth ()
{
    Computation::computeTimeStepWidth ();
    
    double dtAll;
    MPI_Allreduce(&dt_, &dtAll, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    dt_ = dtAll;
};

void ComputationParallel::computePreliminaryVelocities ()
{
    Computation::computePreliminaryVelocities ();
    velocity_communication(discretization_->f, discretization_->g);
};

void ComputationParallel::computeVelocities()
{
    Computation::computeVelocities();
    velocity_communication(discretization_->u, discretization_->v);
}

void ComputationParallel::velocity_communication(double (*u_or_f)(int, int), double (*v_or_g)(int, int))
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
    
    // Communication horizontal first, then vertical!!
    
    // HORIZONTAL COMMUNICATION
    std::vector<MPI_Request> requests_horizontal;
    
    // F Communication (Send)
    // communicate to right
    discretization_->send_boundary_horizontal(to_right, discretization_->uIEnd()-2, discretization_->uJBegin(), discretization_->uJEnd(), discretization_->rank_neighbor(right), u_or_f, discretization_->is_boundary(right));
    // communicate to left
    discretization_->send_boundary_horizontal(to_left, discretization_->uIBegin(), discretization_->uJBegin(), discretization_->uJEnd(), discretization_->rank_neighbor(left), u_or_f, discretization_->is_boundary(left));
    
    // G Communication (Send)
    // communicate to right
    discretization_->send_boundary_horizontal(to_right + 4, discretization_->vIEnd()-1, discretization_->vJBegin(), discretization_->vJEnd()-1, discretization_->rank_neighbor(right), v_or_g, discretization_->is_boundary(right));
    // communicate to left
    discretization_->send_boundary_horizontal(to_left + 4, discretization_->vIBegin(), discretization_->vJBegin(), discretization_->vJEnd()-1, discretization_->rank_neighbor(left), v_or_g, discretization_->is_boundary(left));

    // Receiving F
    MPI_Request current_request;
    // from right
    current_request = discretization_->receive_boundary_horizontal(to_left, discretization_->uIEnd()-1, discretization_->uJBegin(), discretization_->uJEnd(), discretization_->rank_neighbor(right), u_or_f, discretization_->is_boundary(right));
    requests_horizontal.push_back(current_request);
    // from left
    current_request = discretization_->send_boundary_horizontal(to_right, discretization_->uIBegin()-1, discretization_->uJBegin(), discretization_->uJEnd(), discretization_->rank_neighbor(left), u_or_f, discretization_->is_boundary(left)););
    requests_horizontal.push_back(current_request);
    
    // Receiving G
    MPI_Request current_request;
    // from right
    current_request = discretization_->receive_boundary_horizontal(to_left, discretization_->vIEnd(), discretization_->vJBegin(), discretization_->vJEnd()-1, discretization_->rank_neighbor(right), v_or_g, discretization_->is_boundary(right));
    requests_horizontal.push_back(current_request);
    // from left
    current_request = discretization_->send_boundary_horizontal(to_right, discretization_->vIBegin()-1, discretization_->vJBegin(), discretization_->vJEnd(), discretization_->rank_neighbor(left), v_or_g, discretization_->is_boundary(left)););
    requests_horizontal.push_back(current_request);
    
    MPI_Wait(&requests_horizontal,MPI_STATUS_IGNORE)
    
    
    // VERTICAL COMMUNICATION
    std::vector<MPI_Request> requests_vertical;
    
    // F Communication (Send)
    // communicate to below
    discretization_->send_boundary_vertical(down, discretization_->uJBegin(), discretization_->uIBegin(), discretization_->uIEnd(), discretization_->rank_neighbor(below), u_or_f, discretization_->is_boundary(below));
    // communicate to above
    discretization_->send_boundary_vertical(up, discretization_->uJEnd()-1, discretization_->uIBegin(), discretization_->uIEnd(), discretization_->rank_neighbor(above), u_or_f, discretization_->is_boundary(above));
    
    // G Communication (send)
    // communicate to below
    discretization_->send_boundary_vertical(down + 4, discretization_->vJBegin(), discretization_->vIBegin()-1, discretization_->vIEnd(), discretization_->rank_neighbor(below), v_or_g, discretization_->is_boundary(below));
    // communicate to above
    discretization_->send_boundary_vertical(up + 4, discretization_->vJEnd()-2, discretization_->vIBegin()-1, discretization_->vIEnd(), discretization_->rank_neighbor(above), v_or_g, discretization_->is_boundary(above));

    // F Communication (receive)
    // from below
    current_request = discretization_->receive_boundary_vertical(up, discretization_->uJBegin()-1, discretization_->uIBegin(), discretization_->uIEnd(), discretization_->rank_neighbor(below), u_or_f, discretization_->is_boundary(below));
    requests_vertical.push_back(current_request);
    // from above
    current_request = discretization_->receive_boundary_vertical(down, discretization_->uJEnd(), discretization_->uIBegin(), discretization_->uIEnd(), discretization_->rank_neighbor(above), u_or_f, discretization_->is_boundary(above));
    requests_vertical.push_back(current_request);
    
    // G Communication (receive)
    // from below
    current_request = discretization_->receive_boundary_vertical(up, discretization_->vJBegin()-1, discretization_->vIBegin()-1, discretization_->vIEnd(), discretization_->rank_neighbor(below), v_or_g, discretization_->is_boundary(below));
    requests_vertical.push_back(current_request);
    // from above
    current_request = discretization_->receive_boundary_vertical(down, discretization_->vJEnd()-1, discretization_->vIBegin()-1, discretization_->vIEnd(), discretization_->rank_neighbor(above), v_or_g, discretization_->is_boundary(above));
    requests_vertical.push_back(current_request);

    
    MPI_Wait(&requests_vertical,MPI_STATUS_IGNORE)
};    

