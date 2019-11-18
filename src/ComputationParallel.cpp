
#include "ComputationParallel.h"
#include "Partitioning.h"

#include <math.h>
#include <cassert>
#include <assert.h>
#include <vector>

#include <iostream>

void ComputationParallel::initialize (int argc, char *argv[])
{
  //Computation::initialize (int argc, char *argv[]);
  settings_.loadFromFile(argv[1]);
  std::array<int,2> nCellsGlobal = settings_.nCells;
	// computing meshWidth (everywhere with global values)
	double dx = settings_.physicalSize[0]/nCellsGlobal[0];
	double dy = settings_.physicalSize[1]/nCellsGlobal[1];
	meshWidth_ = {dx, dy};
  Partitioning parti;

    // Number of processes determined by command line
    // Get the number of processes
    int MPI_n_processes = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &MPI_n_processes);
    // Get the rank of the process
    int MPI_rank =0;
    MPI_Comm_rank(MPI_COMM_WORLD, &MPI_rank);
    std::cout << "here we are" << MPI_rank << std::endl;

    // compute nCells and assign physical relationships
    if(MPI_rank == 0)
    {
        std::cout << "rank" << MPI_rank << std::endl;
        settings_.printSettings();

        // Computing numer of partitions in each direction
        double ratio_nCellsGlobal = nCellsGlobal[0]/nCellsGlobal[1]; // x/y
        double ratio;
        double ratio_difference;
        double ratio_difference_best = 1000;
        int divisor_best = 1;

        for (int d = 1; d <= MPI_n_processes; d++)
        {
            if (MPI_n_processes % d == 0)
            {
                ratio = d*d/MPI_n_processes;
                ratio_difference = abs(ratio - ratio_nCellsGlobal);
                if (ratio_difference < ratio_difference_best)
                {
                    ratio_difference_best = ratio_difference;
                    divisor_best = d;
                }
            }
        }
        // Number of partitions in both dimensions
        int n_pars_x = divisor_best;
        int n_pars_y = MPI_n_processes/divisor_best;

        // Number of cells in each partiotion except the last one
        int n_Cells_sub_x = int (nCellsGlobal[0] / n_pars_x);
        int n_Cells_sub_y = int (nCellsGlobal[1] / n_pars_y);
        // Number of cells in the last partitions in each dimension
        int n_Cells_sub_x_last = nCellsGlobal[0] - (n_pars_x-1)*n_Cells_sub_x;
        int n_Cells_sub_y_last = nCellsGlobal[1] - (n_pars_y-1)*n_Cells_sub_y;
        // Defining ranks for all partitions
        int ranks_domain[n_pars_x][n_pars_y];
        int rank_counter = 0;
        for(int j = 0; j < n_pars_y; j++)
        {
            for (int i = 0; i < n_pars_x; i++)
            {
                ranks_domain[i][j] = rank_counter;
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
                    std::array<int,4> ranks_neighbors = {ranks_domain[i][j-1], ranks_domain[i+1][j], ranks_domain[i][j+1], ranks_domain[i-1][j]}; // bottom, right, upper, left; caution: check for limits (boundaries)!
                    std::array<bool,4> is_boundary ={ (j-1) == -1, (i+1) == n_pars_x, (j+1) == n_pars_y, (i-1) == -1};
                    std::array<int,2> nCells_sub = {n_Cells_sub_x,n_Cells_sub_y};
                    std::array<int,2> nodeOffset = {n_Cells_sub_x*i, n_Cells_sub_y*j};

                    // setting nCellsGlobal and overwriting nCells for rank 0
                    parti = Partitioning(MPI_rank, ranks_neighbors, is_boundary, nCells_sub, nCellsGlobal, nodeOffset);
                    // discretization_->set_partitioning(MPI_rank, ranks_neighbors, is_boundary, nCells_sub, nCellsGlobal);
                    // Not changing physicalSize because it's not used. Caution: Inconsistent to nCells
                    settings_.nCells = nCells_sub;
                }
                else
                {
                    std::array<int,4> ranks_neighbors = {ranks_domain[i][j-1], ranks_domain[i+1][j], ranks_domain[i][j+1], ranks_domain[i-1][j]}; // bottom, right, upper, left; caution: check for limits (boundaries)!
                    std::array<bool,4> is_boundary = { (j-1) == -1, (i+1) == n_pars_x, (j+1) == n_pars_y, (i-1) == -1};
                    std::array<int,2> nCells_sub = {0,0};
                    std::array<int,2> nodeOffset = {n_Cells_sub_x*i, n_Cells_sub_y*j};
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
                    MPI_Request current_request;
                    MPI_Isend(&ranks_neighbors, 4, MPI_INT, ranks_domain[i][j],0,MPI_COMM_WORLD, &current_request);
                    MPI_Request_free(&current_request);
                    MPI_Isend(&is_boundary, 4, MPI_INT, ranks_domain[i][j],1,MPI_COMM_WORLD, &current_request);
                    MPI_Request_free(&current_request);
                    MPI_Isend(&nCells_sub, 2, MPI_INT, ranks_domain[i][j],2,MPI_COMM_WORLD, &current_request);
                    MPI_Request_free(&current_request);
                   // hier sent einf�gen f�r nodeOffset
                    MPI_Isend(&nodeOffset, 2, MPI_INT, ranks_domain[i][j], 3, MPI_COMM_WORLD, &current_request);
                    MPI_Request_free(&current_request);
                }
            }
        }

    }
    else
    {
        std::array<int,4> ranks_neighbors; // bottom, right, upper, left; caution: check for limits (boundaries)!
        std::array<bool,4> is_boundary;
        std::array<int,2> nCells_sub;
        std::array<int,2> nodeOffset;
        std::vector<MPI_Request> requests;
        MPI_Request current_request;
        MPI_Irecv(&ranks_neighbors, 4, MPI_INT, 0, 0, MPI_COMM_WORLD, &current_request);
        requests.push_back(current_request);
        MPI_Irecv(&is_boundary, 4, MPI_INT, 0, 1, MPI_COMM_WORLD, &current_request);
        requests.push_back(current_request);
        MPI_Irecv(&nCells_sub, 2, MPI_INT, 0, 2, MPI_COMM_WORLD, &current_request);
        requests.push_back(current_request);
        MPI_Irecv(&nodeOffset, 2, MPI_INT, 0, 3, MPI_COMM_WORLD, &current_request);
        requests.push_back(current_request);
        MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);

        // setting nCellsGlobal and overwriting nCells for all ranks except 0
        parti = Partitioning(MPI_rank, ranks_neighbors, is_boundary, nCells_sub, nCellsGlobal, nodeOffset);
        // discretization_->set_partitioning(MPI_rank, ranks_neighbors, is_boundary, nCells_sub, nCellsGlobal);
        // Not changing physicalSize because it's not used. Caution: Inconsistent to nCells
        settings_.nCells = nCells_sub;
    }

    //select DonorCell or CentralDifferences
  	if (settings_.useDonorCell == true)
  	{
  		discretization_ = std::make_shared<DonorCell>(settings_.nCells, meshWidth_, settings_.alpha, parti);
  	}
  	else
  	{
  		discretization_ = std::make_shared<CentralDifferences>(settings_.nCells, meshWidth_,  parti);
  	}

  	//select SOR or GaussSeidel
  	if (settings_.pressureSolver == "SOR")
  	{
  		pressureSolver_ = std::make_unique<SORRedBlack>(discretization_, settings_.epsilon,
  		 settings_.maximumNumberOfIterations, settings_.omega);
  	}
  	else
  	{
  		pressureSolver_ = std::make_unique<GaussSeidel>(discretization_, settings_.epsilon,
  		 settings_.maximumNumberOfIterations);
  	}

  	//initialize outputWriters
   	outputWriterParaviewParallel_ = std::make_unique<OutputWriterParaviewParallel>(discretization_, parti);
  	outputWriterTextParallel_ = std::make_unique<OutputWriterTextParallel>(discretization_, parti);

    std::cout << "rank" << MPI_rank << std::endl;
    settings_.printSettings();


};

void ComputationParallel::runSimulation ()
{
	double time = 0;
  int MPI_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &MPI_rank);
  std::cout << "here we are 2: " << MPI_rank << std::endl;
	while(time<settings_.endTime)
	{
		applyBoundaryValues();
    std::cout << "here we are 3: " << MPI_rank << std::endl;

		// std::cout << "time" << time << std::endl;
		// compute dt_ and time
		computeTimeStepWidth();
    std::cout << "here we are 4: " << MPI_rank << std::endl;

		if(time+dt_>settings_.endTime) dt_ = settings_.endTime - time;
		time += dt_;
		//std::cout << "time_step" << dt_ << std::endl;

    std::cout << "here we are 5: " << MPI_rank << std::endl;


		// compute f and g
		computePreliminaryVelocities();
		// outputWriterText_->writeFile(time);
    std::cout << "here we are 6: " << MPI_rank << std::endl;

		//compute rhs
		computeRightHandSide();
		//compute p with SOR or GaussSeidel

    std::cout << "here we are 7: " << MPI_rank << std::endl;


		computePressure();
		//compute u and v

    std::cout << "here we are 8: " << MPI_rank << std::endl;

		computeVelocities();

    std::cout << "here we are 9: " << MPI_rank << std::endl;


		outputWriterParaviewParallel_->writeFile(time);
		outputWriterTextParallel_->writeFile(time);
		outputWriterTextParallel_->writePressureFile();
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
  int MPI_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &MPI_rank);
  std::cout << "here we are 5.1: " << MPI_rank << std::endl;
    Computation::computePreliminaryVelocities ();
      std::cout << "here we are 5.2: " << MPI_rank << std::endl;
    preliminaryVelocity_communication();
      std::cout << "here we are 5.3: " << MPI_rank << std::endl;
};

void ComputationParallel::computeVelocities()
{
    Computation::computeVelocities();
    finalVelocity_communication();
}

void ComputationParallel::finalVelocity_communication() //NUR Kommunikation von f und g, was ist mit u und v???
{
  int MPI_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &MPI_rank);
  std::cout << "here we are 6.2.1: " << MPI_rank << std::endl;
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

    // Communication horizontal first, then vertical!! // warum?

    // HORIZONTAL COMMUNICATION
    std::vector<MPI_Request> requests_horizontal;

    // F Communication (Send)
    // communicate to right
    discretization_->send_boundary_horizontal_u(to_right, discretization_->uIEnd()-2, discretization_->uJBegin(), discretization_->uJEnd(), discretization_->rank_neighbor(right), discretization_->is_boundary(right));
    // communicate to left
    discretization_->send_boundary_horizontal_u(to_left, discretization_->uIBegin(), discretization_->uJBegin(), discretization_->uJEnd(), discretization_->rank_neighbor(left), discretization_->is_boundary(left));

    // G Communication (Send)
    // communicate to right
    discretization_->send_boundary_horizontal_v(to_right + 4, discretization_->vIEnd()-1, discretization_->vJBegin(), discretization_->vJEnd()-1, discretization_->rank_neighbor(right), discretization_->is_boundary(right));
    // communicate to left
    discretization_->send_boundary_horizontal_v(to_left + 4, discretization_->vIBegin(), discretization_->vJBegin(), discretization_->vJEnd()-1, discretization_->rank_neighbor(left), discretization_->is_boundary(left));

    std::cout << "here we are 6.2.2: " << MPI_rank << std::endl;

    std::cout << "request_length: " << requests_horizontal.size() << MPI_rank << std::endl;


    // Receiving F
    MPI_Request current_request;
    // from right
    discretization_->receive_boundary_horizontal_u(current_request, to_left, discretization_->uIEnd()-1, discretization_->uJBegin(), discretization_->uJEnd(), discretization_->rank_neighbor(right), discretization_->is_boundary(right));
    std::cout << "is_boundary(right): " << discretization_->is_boundary(right) << MPI_rank << std::endl;
    if(!discretization_->is_boundary(right)) requests_horizontal.push_back(current_request);

    std::cout << "request_length: " << requests_horizontal.size() << MPI_rank << std::endl;

    // from left
    discretization_->receive_boundary_horizontal_u(current_request, to_right, discretization_->uIBegin()-1, discretization_->uJBegin(), discretization_->uJEnd(), discretization_->rank_neighbor(left), discretization_->is_boundary(left));
    if(!discretization_->is_boundary(left)) requests_horizontal.push_back(current_request);

    std::cout << "request_length: " << requests_horizontal.size() << MPI_rank << std::endl;


    // Receiving G
    // from right
    discretization_->receive_boundary_horizontal_v(current_request, to_left + 4, discretization_->vIEnd(), discretization_->vJBegin(), discretization_->vJEnd()-1, discretization_->rank_neighbor(right), discretization_->is_boundary(right));
    if(!discretization_->is_boundary(right)) requests_horizontal.push_back(current_request);
    // from left
    discretization_->receive_boundary_horizontal_v(current_request, to_right + 4, discretization_->vIBegin()-1, discretization_->vJBegin(), discretization_->vJEnd(), discretization_->rank_neighbor(left), discretization_->is_boundary(left));
    if(!discretization_->is_boundary(left)) requests_horizontal.push_back(current_request);

    std::cout << "here we are 6.2.3: " << MPI_rank << std::endl;

    std::cout << "request_length: " << requests_horizontal.size() << MPI_rank << std::endl;

    MPI_Waitall(requests_horizontal.size(), requests_horizontal.data() ,MPI_STATUSES_IGNORE);

    std::cout << "here we are 6.2.4: " << MPI_rank << std::endl;

    // VERTICAL COMMUNICATION
    std::vector<MPI_Request> requests_vertical;

    // F Communication (Send)
    // communicate to below
    discretization_->send_boundary_vertical_u(down, discretization_->uJBegin(), discretization_->uIBegin(), discretization_->uIEnd(), discretization_->rank_neighbor(below), discretization_->is_boundary(below));
    // communicate to above
    discretization_->send_boundary_vertical_u(up, discretization_->uJEnd()-1, discretization_->uIBegin(), discretization_->uIEnd(), discretization_->rank_neighbor(above), discretization_->is_boundary(above));

    // G Communication (send)
    // communicate to below
    discretization_->send_boundary_vertical_g(down + 4, discretization_->vJBegin(), discretization_->vIBegin()-1, discretization_->vIEnd(), discretization_->rank_neighbor(below), discretization_->is_boundary(below));
    // communicate to above
    discretization_->send_boundary_vertical_g(up + 4, discretization_->vJEnd()-2, discretization_->vIBegin()-1, discretization_->vIEnd(), discretization_->rank_neighbor(above), discretization_->is_boundary(above));

    // F Communication (receive)
    // from below
    discretization_->receive_boundary_vertical_u(current_request, up, discretization_->uJBegin()-1, discretization_->uIBegin(), discretization_->uIEnd(), discretization_->rank_neighbor(below), discretization_->is_boundary(below));
    if(!discretization_->is_boundary(below)) requests_vertical.push_back(current_request);
    // from above
    discretization_->receive_boundary_vertical_u(current_request, down, discretization_->uJEnd(), discretization_->uIBegin(), discretization_->uIEnd(), discretization_->rank_neighbor(above), discretization_->is_boundary(above));
    if(!discretization_->is_boundary(above)) requests_vertical.push_back(current_request);

    // G Communication (receive)
    // from below
    discretization_->receive_boundary_vertical_g(current_request, up + 4, discretization_->vJBegin()-1, discretization_->vIBegin()-1, discretization_->vIEnd(), discretization_->rank_neighbor(below), discretization_->is_boundary(below));
    if(!discretization_->is_boundary(below)) requests_vertical.push_back(current_request);
    // from above
    discretization_->receive_boundary_vertical_g(current_request, down + 4, discretization_->vJEnd()-1, discretization_->vIBegin()-1, discretization_->vIEnd(), discretization_->rank_neighbor(above), discretization_->is_boundary(above));
    if(!discretization_->is_boundary(above)) requests_vertical.push_back(current_request);
    std::cout << "here we are 6.2.5: " << MPI_rank << std::endl;


    MPI_Waitall(requests_vertical.size(), requests_vertical.data(),MPI_STATUSES_IGNORE);

    std::cout << "here we are 6.2.6: " << MPI_rank << std::endl;

};

void ComputationParallel::preliminaryVelocity_communication() //NUR Kommunikation von f und g, was ist mit u und v???
{
  int MPI_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &MPI_rank);
  std::cout << "here we are 5.2.1: " << MPI_rank << std::endl;
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

    // Communication horizontal first, then vertical!! // warum?

    // HORIZONTAL COMMUNICATION
    std::vector<MPI_Request> requests_horizontal;

    // F Communication (Send)
    // communicate to right
    discretization_->send_boundary_horizontal_f(to_right, discretization_->uIEnd()-2, discretization_->uJBegin(), discretization_->uJEnd(), discretization_->rank_neighbor(right), discretization_->is_boundary(right));
    // communicate to left
    discretization_->send_boundary_horizontal_f(to_left, discretization_->uIBegin(), discretization_->uJBegin(), discretization_->uJEnd(), discretization_->rank_neighbor(left), discretization_->is_boundary(left));

    // G Communication (Send)
    // communicate to right
    discretization_->send_boundary_horizontal_g(to_right + 4, discretization_->vIEnd()-1, discretization_->vJBegin(), discretization_->vJEnd()-1, discretization_->rank_neighbor(right), discretization_->is_boundary(right));
    // communicate to left
    discretization_->send_boundary_horizontal_g(to_left + 4, discretization_->vIBegin(), discretization_->vJBegin(), discretization_->vJEnd()-1, discretization_->rank_neighbor(left), discretization_->is_boundary(left));

    std::cout << "here we are 5.2.2: " << MPI_rank << std::endl;

    std::cout << "request_length: " << requests_horizontal.size() << MPI_rank << std::endl;


    // Receiving F
    MPI_Request current_request;
    // from right
    discretization_->receive_boundary_horizontal_f(current_request, to_left, discretization_->uIEnd()-1, discretization_->uJBegin(), discretization_->uJEnd(), discretization_->rank_neighbor(right), discretization_->is_boundary(right));
    std::cout << "is_boundary(right): " << discretization_->is_boundary(right) << MPI_rank << std::endl;
    if(!discretization_->is_boundary(right)) requests_horizontal.push_back(current_request);

    std::cout << "request_length: " << requests_horizontal.size() << MPI_rank << std::endl;

    // from left
    discretization_->receive_boundary_horizontal_f(current_request, to_right, discretization_->uIBegin()-1, discretization_->uJBegin(), discretization_->uJEnd(), discretization_->rank_neighbor(left), discretization_->is_boundary(left));
    if(!discretization_->is_boundary(left)) requests_horizontal.push_back(current_request);

    std::cout << "request_length: " << requests_horizontal.size() << MPI_rank << std::endl;


    // Receiving G
    // from right
    discretization_->receive_boundary_horizontal_g(current_request, to_left + 4, discretization_->vIEnd(), discretization_->vJBegin(), discretization_->vJEnd()-1, discretization_->rank_neighbor(right), discretization_->is_boundary(right));
    if(!discretization_->is_boundary(right)) requests_horizontal.push_back(current_request);
    // from left
    discretization_->receive_boundary_horizontal_g(current_request, to_right + 4, discretization_->vIBegin()-1, discretization_->vJBegin(), discretization_->vJEnd(), discretization_->rank_neighbor(left), discretization_->is_boundary(left));
    if(!discretization_->is_boundary(left)) requests_horizontal.push_back(current_request);

    std::cout << "here we are 5.2.3: " << MPI_rank << std::endl;

    std::cout << "request_length: " << requests_horizontal.size() << MPI_rank << std::endl;

    MPI_Waitall(requests_horizontal.size(), requests_horizontal.data() ,MPI_STATUSES_IGNORE);

    std::cout << "here we are 5.2.4: " << MPI_rank << std::endl;

    // VERTICAL COMMUNICATION
    std::vector<MPI_Request> requests_vertical;

    // F Communication (Send)
    // communicate to below
    discretization_->send_boundary_vertical_f(down, discretization_->uJBegin(), discretization_->uIBegin(), discretization_->uIEnd(), discretization_->rank_neighbor(below), discretization_->is_boundary(below));
    // communicate to above
    discretization_->send_boundary_vertical_f(up, discretization_->uJEnd()-1, discretization_->uIBegin(), discretization_->uIEnd(), discretization_->rank_neighbor(above), discretization_->is_boundary(above));

    // G Communication (send)
    // communicate to below
    discretization_->send_boundary_vertical_g(down + 4, discretization_->vJBegin(), discretization_->vIBegin()-1, discretization_->vIEnd(), discretization_->rank_neighbor(below), discretization_->is_boundary(below));
    // communicate to above
    discretization_->send_boundary_vertical_g(up + 4, discretization_->vJEnd()-2, discretization_->vIBegin()-1, discretization_->vIEnd(), discretization_->rank_neighbor(above), discretization_->is_boundary(above));

    // F Communication (receive)
    // from below
    discretization_->receive_boundary_vertical_f(current_request, up, discretization_->uJBegin()-1, discretization_->uIBegin(), discretization_->uIEnd(), discretization_->rank_neighbor(below), discretization_->is_boundary(below));
    if(!discretization_->is_boundary(below)) requests_vertical.push_back(current_request);
    // from above
    discretization_->receive_boundary_vertical_f(current_request, down, discretization_->uJEnd(), discretization_->uIBegin(), discretization_->uIEnd(), discretization_->rank_neighbor(above), discretization_->is_boundary(above));
    if(!discretization_->is_boundary(above)) requests_vertical.push_back(current_request);

    // G Communication (receive)
    // from below
    discretization_->receive_boundary_vertical_g(current_request, up + 4, discretization_->vJBegin()-1, discretization_->vIBegin()-1, discretization_->vIEnd(), discretization_->rank_neighbor(below), discretization_->is_boundary(below));
    if(!discretization_->is_boundary(below)) requests_vertical.push_back(current_request);
    // from above
    discretization_->receive_boundary_vertical_g(current_request, down + 4, discretization_->vJEnd()-1, discretization_->vIBegin()-1, discretization_->vIEnd(), discretization_->rank_neighbor(above), discretization_->is_boundary(above));
    if(!discretization_->is_boundary(above)) requests_vertical.push_back(current_request);
    std::cout << "here we are 5.2.5: " << MPI_rank << std::endl;


    MPI_Waitall(requests_vertical.size(), requests_vertical.data(),MPI_STATUSES_IGNORE);

    std::cout << "here we are 5.2.6: " << MPI_rank << std::endl;

};



void ComputationParallel::applyBoundaryValues ()
{
  // neigbor indices
  int below = 0;
  int above = 2;
  int right = 1;
  int left = 3;

	// u,f setting
	// lower u ghost layer without corners
  if (discretization_->is_boundary(below))
  {
  	int j = discretization_->uJBegin()-1;
  	for(int i = discretization_->uIBegin(); i < discretization_->uIEnd()-1; i++)
  	{
  		discretization_->u(i,j) = 2*settings_.dirichletBcBottom[0]-discretization_->u(i,j+1);
  		discretization_->f(i,j) = discretization_->u(i,j);
  	};
  }
	// upper u ghost layer without corners
  if (discretization_->is_boundary(above))
  {
  	int j = discretization_->uJEnd();
  	for(int i = discretization_->uIBegin(); i < discretization_->uIEnd()-1; i++)
  	{
  		discretization_->u(i,j) = 2*settings_.dirichletBcTop[0]-discretization_->u(i,j-1);
  		discretization_->f(i,j) = discretization_->u(i,j);
  	};
  }
	// left u ghost layer with corners
  if (discretization_->is_boundary(left))
  {
  	int i = discretization_->uIBegin()-1;
  	for(int j = discretization_->uJBegin()-1; j < discretization_->uJEnd()+1; j++)
  	{
  		discretization_->u(i,j) = settings_.dirichletBcLeft[0];
  		discretization_->f(i,j) = discretization_->u(i,j);
  	}
  }
	// right u Nathi-not ghost layer with corners
  if (discretization_->is_boundary(right))
  {
  	int i = discretization_->uIEnd()-1;
  	for(int j = discretization_->uJBegin()-1; j < discretization_->uJEnd()+1; j++)
  	{
  		discretization_->u(i,j) = settings_.dirichletBcRight[0];
  		discretization_->f(i,j) = discretization_->u(i,j);
  	}
  }

	// v,g setting
	// lower v ghost layer without corners
  if (discretization_->is_boundary(below))
  {
  	int j = discretization_->vJBegin()-1;
  	for(int i = discretization_->vIBegin(); i < discretization_->vIEnd(); i++)
  	{
  		discretization_->v(i,j) = settings_.dirichletBcBottom[1];
  		discretization_->g(i,j) = discretization_->v(i,j);
  	};
  }
	// upper v  Nathi-not ghost layer without corners
  if (discretization_->is_boundary(above))
  {
  	int j = discretization_->vJEnd()-1;
  	for(int i = discretization_->vIBegin(); i < discretization_->vIEnd(); i++)
  	{
  		discretization_->v(i,j) = settings_.dirichletBcTop[1];
  		discretization_->g(i,j) = discretization_->v(i,j);
  	};
  }
	// left v ghost layer with corners
  if (discretization_->is_boundary(left))
  {
  	int i = discretization_->vIBegin()-1;
  	for(int j = discretization_->vJBegin()-1; j < discretization_->vJEnd(); j++)
  	{
  		discretization_->v(i,j) = 2*settings_.dirichletBcLeft[1]-discretization_->v(i+1,j);
  		discretization_->g(i,j) = discretization_->v(i,j);
  	}
  }
	// right v ghost layer with corners
  if (discretization_->is_boundary(right))
  {
  	int i = discretization_->vIEnd();
  	for(int j = discretization_->vJBegin()-1; j < discretization_->vJEnd(); j++)
  	{
  		discretization_->v(i,j) = 2*settings_.dirichletBcRight[1]-discretization_->v(i-1,j);
  		discretization_->g(i,j) = discretization_->v(i,j);
  	}
  }
};

//void ComputationParallel::computeRightHandSide ????
