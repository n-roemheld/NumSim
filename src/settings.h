
#pragma once

#include <iostream>
#include <array>
#include "PressureSolver.h"
#include "Discretization.h"
#include "Cycle.h"
#include <fstream>   // for file operations
#include <iostream>  // for cout

/** All settings that parametrize a simulation run.
 */
struct Settings
{
  std::array<int,2> nCells;          //< number of cells in x and y direction
  std::array<double,2> physicalSize; //< physical size of the domain
  double re = 1000;                  //< reynolds number
  double endTime = 10.0;             //< end time of the simulation
  double tau = 0.5;                  //< safety factor for time step width
  double maximumDt = 0.1;            //< maximum time step width

  std::array<double,2> g;    //< external forces

  bool useDonorCell = false;         //< if the donor cell scheme schould be used
  double alpha = 0.5;                //< factor for donor-cell scheme

  std::array<double,2> dirichletBcBottom;  //< prescribed values of u,v at bottom of domain
  std::array<double,2> dirichletBcTop;     //< prescribed values of u,v at top of domain
  std::array<double,2> dirichletBcLeft;    //< prescribed values of u,v at left of domain
  std::array<double,2> dirichletBcRight;   //< prescribed values of u,v at right of domain

  std::string pressureSolver = "SOR";      //< which pressure solver to use, "GaussSeidel" or "SOR" or "multigrid"
  double omega = 1.0;                //< overrelaxation factor
  double epsilon = 1e-5;             //< tolerance for the residual in the pressure solver
  int maximumNumberOfIterations = 1e5;    //< maximum number of iterations in the solver

  // multigrid parameters
  std::string smoother;
  std::string coarser;
  std::string endSolver;
  std::string cycle;

  int numberOfIterationsPre;
  int numberOfIterationsPost;

  bool recursive;
  int maxLevel;
  std::vector<int> gamma; // todo: add to load from file

  //! parse a text file with settings, each line contains "<parameterName> = <value>"
  void loadFromFile(std::string filename);

  //! output all settings to console
  void printSettings();
};
