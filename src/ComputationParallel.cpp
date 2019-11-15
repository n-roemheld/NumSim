
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
                    std::vector<int,4> ranks_neighbors = {ranks_domain(i,j-1), ranks_domain(i+1,j), ranks_domain(i,j+1), ranks_domain(i-1,j)}; // bottom, right, upper, left; caution: check for limits (boundaries)!
                    std::vector<int,4> is_boundary = { (j-1) == -1, (i+1) == n_pars_x, (j+1) == n_pars_y, (i-1) == -1};
                    std::vector<int,2> nCells_sub = {n_Cells_sub_x,n_Cells_sub_y};
                    discretization_->set_partitioning(MPI_rank, ranks_neighbors, is_boundary, nCells);
                }
                else
                {
                    std::vector<int,4> ranks_neighbors = {ranks_domain(i,j-1), ranks_domain(i+1,j), ranks_domain(i,j+1), ranks_domain(i-1,j)}; // bottom, right, upper, left; caution: check for limits (boundaries)!
                    std::vector<int,4> is_boundary = { (j-1) == -1, (i+1) == n_pars_x, (j+1) == n_pars_y, (i-1) == -1};
                    std::vector<int,2> nCells_sub = {0,0};
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
        std::vector<int,4> ranks_neighbors; // bottom, right, upper, left; caution: check for limits (boundaries)!
        std::vector<int,4> is_boundary;
        std::vector<int,2> nCells_sub;
        MPI_Irecv(&ranks_neighbors, 10, MPI_INT, 0, 0, MPI_COMM_WORLD);
        MPI_Irecv(&is_boundary, 10, MPI_INT, 0, 1, MPI_COMM_WORLD);
        MPI_Irecv(&nCells_sub, 10, MPI_INT, 0, 2, MPI_COMM_WORLD);
        discretization_->set_partitioning(MPI_rank, ranks_neighbors, is_boundary, nCells_sub);
    }
    MPI_Wait()
};

void ComputationParallel::computeTimeStepWidth ()
{
    Computation::computeTimeStepWidth ();
    
    
    MPI_Allreduce(&dt_, &dtAll_, partitioning_.getSize(), MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    dt_ = dtAll_;
};