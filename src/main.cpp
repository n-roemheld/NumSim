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
  ComputationParallel comp;
  comp.initialize(argc, argv);
  time_t startTime;
  time_t endTime;
  std::cout << "StartTime: " << time(&startTime) << std::endl;
  comp.runSimulation();
  std::cout << "EndTime: "<< time(&endTime) << std::endl;
  std::cout << "Difference: " << -(startTime-endTime) << std::endl;

  std::cout << "here we are 130: "  << std::endl;


  MPI_Finalize();
  std::cout << "here we are 131: "  << std::endl;

  return EXIT_SUCCESS;
}
