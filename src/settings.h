
#pragma once

#include "FieldVariable.h"
#include <iostream>
#include <array>
#include <memory>

/** All settings that parametrize a simulation run.
 */
struct Settings
{
  std::array<int,2> nCells;          //< number of cells in x and y direction
  std::array<double,2> physicalSize; //< physical size of the domain
  double re = 1000;                  //< reynolds number
  double prandtl = 0.71;             //< prandtl number
  double beta = 1.019367991845056;   //< Volumetric expansion coefficient beta = Ri / (g * Re^2)
  double endTime = 10.0;             //< end time of the simulation
  double tau = 0.5;                  //< safety factor for time step width
  double maximumDt = 0.1;            //< maximum time step width

  double outputFileEveryDt = 0;        //< call ouput writers only every dt (or at every iteration if provided value is too small)

  std::array<double,2> g;    //< external forces

  std::string participantName = "NumSim";
  std::string meshName = "NumSimMesh";
  std::string preciceConfigFile = "precice-config.xml";
  std::string readDataName = "Heat-Flux";
  std::string writeDataName = "Temperature";

  int vertexSize = 0;
	std::vector<int> vertex_i;
	std::vector<int> vertex_j;
  std::vector<double> coords;

  bool useDonorCell = false;         //< if the donor cell scheme schould be used
  double alpha = 0.5;                //< factor for donor-cell scheme for pressure
  double gamma = 0.5;				 //< factor for donor-cell scheme for temperature

  double xOrigin = 0; // Origin for the geometry (bottom left corner)
  double yOrigin = 0;

  std::array<double,2> dirichletBcBottom;  //< prescribed values of u,v at bottom of domain
  std::array<double,2> dirichletBcTop;     //< prescribed values of u,v at top of domain
  std::array<double,2> dirichletBcLeft;    //< prescribed values of u,v at left of domain
  std::array<double,2> dirichletBcRight;   //< prescribed values of u,v at right of domain

  std::string pressureSolver = "SOR";      //< which pressure solver to use, "GaussSeidel" or "SOR"
  double omega = 1.0;                //< overrelaxation factor
  double epsilon = 1e-5;             //< tolerance for the residual in the pressure solver
  int maximumNumberOfIterations = 1e5;    //< maximum number of iterations in the solver

  double uInit_; //< initial value for velocity u
  double vInit_; //< initial value for velocity v
  double pInit_; //< initial value for pressure p
  double TInit_; //< initial value for temperature T

  std::shared_ptr<Array2D> geometryPVString_; //< describes typ of cell for pressure
  std::shared_ptr<Array2D> geometryPVOrientation_; //< describes direction of fluid cells; 1,2,3,4 = left,right,lower,upper
  std::shared_ptr<Array2D> geometryPV1_;
  std::shared_ptr<Array2D> geometryPV2_;

  std::shared_ptr<Array2D> geometryTString_;
  std::shared_ptr<Array2D> geometryT1_;

  std::string geometryFile = "";
  std::vector<int> orientation_;


  //! parse a text file with settings, each line contains "<parameterName> = <value>"
  void loadFromFile(std::string filename);

  //! output all settings to console
  void printSettings();

  void loadGeometryFile();
};
