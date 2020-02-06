#include "settings.h"
#include "Computation.h"

#include <iostream>
#include <cstdlib>
#include <vector>
#include <chrono>


int main(int argc, char *argv[])
{
  // if the number of given command line arguments is only 1 (= the program name), print out usage information and exit
  if (argc == 1)
  {
    std::cout << "usage: " << argv[0] << " <filename>" << std::endl;

    return EXIT_FAILURE;
  }

  // creating data independend of MG parameters in input file
  if(argc == 3 && argv[2] == std::string("-data"))
  {
    std::cout << "DATA" << std::endl;
  //------added for creating data----------
  std::vector < std::string > smootherStrings;
  // smootherStrings.push_back("Jacobi");
  // smootherStrings.push_back("DJacobi");
  smootherStrings.push_back("GS");
  std::vector < std::string > coarserStrings;
  coarserStrings.push_back("Default");
  std::vector < std::string > endSolverStrings;
  endSolverStrings.push_back("None");
  // endSolverStrings.push_back("GS");
  std::vector < int > mls;
  for(int i = 3; i <= 3; i++)
  {
    mls.push_back(i);
  }
  std::vector < int > noIPres;
  for(int i = 2; i <= 2; i++)
  {
    noIPres.push_back(i);
  }
  std::vector < int > noIPosts;
  for(int i = 2; i <= 2; i++)
  {
    noIPosts.push_back(i);
  }
  std::vector < std::string > cycleStrings;
  cycleStrings.push_back("V");
  cycleStrings.push_back("W");
  //----------------------------------------


  Computation comp;

  
  //------added for creating data----------
  for(std::string smootherString : smootherStrings)
  {
    for(std::string coarserString : coarserStrings)
    {
      for(std::string endSolverString : endSolverStrings)
      {
        for(int ml : mls)
        {
          for(int noIPre : noIPres)
          {
            for(int noIPost : noIPosts)
            {
              for(std::string cycleString : cycleStrings)
              {
  //---------------------------------------

  //------changed for creating data----------
  comp.initialize(argc, argv, smootherString, coarserString, endSolverString, ml, noIPre, noIPost, cycleString);
  // comp.initialize(argc, argv);
  //----------------------------------------

  time_t startTime;
  time_t endTime;
  // std::cout << "StartTime: " << time(&startTime) << std::endl;
  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
  comp.runSimulation();
  std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
  // std::cout << "EndTime: "<< time(&endTime) << std::endl;
  // std::cout << "Difference: " << -(startTime-endTime) << std::endl;
  std::cout << "Time difference: " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << "[s]" << std::endl;

  
  //------added for creating data----------
              }
            }
          }
        }
      }
    }
  }
  //---------------------------------------
  }
  // run simulation with parameters from input file
  else
  {
    Computation comp;

    comp.initialize(argc, argv, "", "","" , 0, 0, 0, "");

    time_t startTime;
    time_t endTime;
    // std::cout << "StartTime: " << time(&startTime) << std::endl;
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    comp.runSimulation();
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    // std::cout << "EndTime: "<< time(&endTime) << std::endl;
    // std::cout << "Difference: " << -(startTime-endTime) << std::endl;
    std::cout << "Time difference: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;
  }


  return EXIT_SUCCESS;
}
