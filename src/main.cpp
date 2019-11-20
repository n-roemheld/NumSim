#include "settings.h"
#include "ComputationParallel.h"

#include <iostream>
#include <cstdlib>

#include <mpi.h>
#include "Partitioning.h"


int main(int argc, char *argv[])
{
  // if the number of given command line arguments is only 1 (= the program name), print out usage information and exit
  if (argc == 1)
  {
    std::cout << "usage: " << argv[0] << " <filename>" << std::endl;

    return EXIT_FAILURE;
  }


  MPI_Init(&argc, &argv);

  int MPI_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &MPI_rank);

  // if(MPI_rank == 1)
  // {
  //   int buffer= 0;
  //   MPI_Request req;
  //   MPI_Isend(&buffer, 1, MPI_INT, 0, 7, MPI_COMM_WORLD, &req);
  // }
  // if(MPI_rank == 0)
  // {
  //   MPI_Request req2;
  //   int buffer2;
  //   MPI_Irecv(&buffer2, 1, MPI_INT, 1, 7, MPI_COMM_WORLD, &req2);
  //
  //   MPI_Wait(&req2, MPI_STATUS_IGNORE);
  // }



  ComputationParallel comp;
  comp.initialize(argc, argv);
  time_t startTime;
  time_t endTime;
  std::cout << "StartTime: " << time(&startTime) << std::endl;
  comp.runSimulation();
  std::cout << "EndTime: "<< time(&endTime) << std::endl;
  std::cout << "Difference: " << -(startTime-endTime) << std::endl;

  std::cout << "here we are 130: "  << MPI_rank << std::endl;


  MPI_Finalize();
  std::cout << "here we are 131: " << MPI_rank  << std::endl;

  return EXIT_SUCCESS;
}
