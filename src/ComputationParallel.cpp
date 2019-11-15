
#include "ComputationParallel.h"
#include "Partitioning.h"

#include <math.h>
#include <cassert>
#include <assert.h>
#include <vector>

void ComputationParallel::initialize (int argc, char *argv[])
{
    Computation::initialize (int argc, char *argv[]);
    
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
                    discretization_->set_partitioning(MPI_rank, ranks_neighbors, is_boundary, nCells);
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
        discretization_->set_partitioning(MPI_rank, ranks_neighbors, is_boundary, nCells_sub);
    }
    MPI_Wait(&requests,MPI_STATUS_IGNORE)
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
    std::vector<MPI_Request> requests;
    
    // F Communication (Send)
    // communicate to below
    int below = 0;
    send_boundary_vertical(below, discretization_->uJBegin(), discretization_->uIBegin(), discretization_->uIEnd()-1, discretization_->rank_neighbor(below), discretization_->f, discretization_->is_boundary(below));
    // communicate to above
    int above = 2;
    send_boundary_vertical(above, discretization_->uJEnd()-1, discretization_->uIBegin(), discretization_->uIEnd()-1, discretization_->rank_neighbor(above), discretization_->f, discretization_->is_boundary(above));
    // communicate to right
    int right = 1;
    send_boundary_horizontal(right, discretization_->uIEnd()-2, discretization_->uJBegin(), discretization_->uJEnd(), discretization_->rank_neighbor(right), discretization_->f, discretization_->is_boundary(right));
    // communicate to left
    int left = 3;
    send_boundary_horizontal(left, discretization_->uIBegin(), discretization_->uJBegin(), discretization_->uJEnd(), discretization_->rank_neighbor(left), discretization_->f, discretization_->is_boundary(left));
    
    // G Communication (Send)
    // communicate to below
    int below = 0;
    send_boundary_vertical(below + 4, discretization_->vJBegin(), discretization_->vIBegin(), discretization_->vIEnd(), discretization_->rank_neighbor(below), discretization_->g, discretization_->is_boundary(below));
    // communicate to above
    int above = 2;
    send_boundary_vertical(above + 4, discretization_->vJEnd()-2, discretization_->vIBegin(), discretization_->vIEnd(), discretization_->rank_neighbor(above), discretization_->g, discretization_->is_boundary(above));
    // communicate to right
    int right = 1;
    send_boundary_horizontal(right + 4, discretization_->vIEnd()-1, discretization_->vJBegin(), discretization_->vJEnd()-1, discretization_->rank_neighbor(right), discretization_->g, discretization_->is_boundary(right));
    // communicate to left
    int left = 3;
    send_boundary_horizontal(left + 4, discretization_->vIBegin(), discretization_->vJBegin(), discretization_->vJEnd()-1, discretization_->rank_neighbor(left), discretization_->g, discretization_->is_boundary(left));
    
    // Receiving F ___________________________________________________________________ ToDo ___________________________________________________________________
    std::vector<MPI_Request> requests;
    MPI_Request current_request;
    // below
    current_request = send_boundary_vertical(int tag, int j_fixed, int i_begin, int i_end, int target_rank, double (*fVar)(int, int), bool do_nothing);
    requests.push_back(current_request);
    // above
    current_request = send_boundary_vertical();
    requests.push_back(current_request);
    // right
    current_request = send_boundary_horizontal();
    requests.push_back(current_request);
    //left
    current_request = send_boundary_horizontal();
    requests.push_back(current_request);
    
    MPI_Wait(&requests,MPI_STATUS_IGNORE)
    
    // for each neighbor
    //      if not_boundary
    //          send_to neighbor
    //          receive_from_neighbor
    //          set_inner_boundaries_according_to_received_data
    //      else
    //          do_nothing!
    //      end
    // end
    // MPI_Wait
};
        

void ComputationParallel::send_boundary_vertical(int tag, int j_fixed, int i_begin, int i_end, int target_rank, double (*fVar)(int, int), bool do_nothing)
{
    if (!do_nothing)
    {
        std::vector<double> send_buffer(i_end - i_begin);
        int bi = 0;
        for(int i = i_begin, i < i_end, i++)
        {
            send_buffer[bi] = fVar(i,j);
            bi++;
        }
        MPI_Isend(&send_buffer, i_end-i_begin, MPI_DOUBLE, target_rank, direction, MPI_COMM_WORLD);

    }
};

void ComputationParallel::send_boundary_horizontal(int tag, int i_fixed, int j_begin, int j_end, int target_rank, double (*fVar)(int, int), bool do_nothing)
{
    if (!do_nothing)
    {
        std::vector<double> send_buffer(j_end - j_begin);
        int bj = 0;
        for(int j = j_begin, j < j_end, j++)
        {
            send_buffer[bj] = fVar(i,j);
            bj++;
        }
        MPI_Isend(&send_buffer, j_end-j_begin, MPI_DOUBLE, target_rank, direction, MPI_COMM_WORLD);

    }
};

MPI_Request ComputationParallel::receive_boundary_vertical(int sender_tag, int j_fixed, int i_begin, int i_end, int source_rank, double (*fVar)(int, int), bool do_nothing)
{
    if (!do_nothing)
    {
        MPI_Request current_request;
        std::vector<double> rcv_buffer;
        MPI_Irecv(&rcv_buffer, i_end - i_begin, MPI_DOUBLE, source_rank, sender_tag, MPI_COMM_WORLD, &current_request);
        int j = j_fixed;
        int bi = 0;
        for (int i = i_begin; i < i_end; i++)
        {
            fVar(i,j) = rcv_buffer(bi);
            bi++;
        }
    }
}

MPI_Request ComputationParallel::receive_boundary_horizontal(int sender_tag, int i_fixed, int j_begin, int j_end, int source_rank, double (*fVar)(int, int), bool do_nothing)
{
    if (!do_nothing)
    {
        MPI_Request current_request;
        std::vector<double> rcv_buffer;
        MPI_Irecv(&rcv_buffer, j_end - j_begin, MPI_DOUBLE, source_rank, sender_tag, MPI_COMM_WORLD, &current_request);
        int i = i_fixed;
        int bj = 0;
        for (int j = j_begin; j < j_end; j++)
        {
            fVar(i,j) = rcv_buffer(bj);
            bj++;
        }
    }
}