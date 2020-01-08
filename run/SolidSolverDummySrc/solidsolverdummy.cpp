/*
 * Copyright (C) 2019  Author: Alexander Jaust
 *
 * Contact: alexander.jaust@ipvs.uni-stuttgart.de
 *
 * No license assigned yet.
 */

#include<algorithm>
#include<array>
#include<cassert>
#include<iostream>
#include "SolverInterface.hpp"

// Domain parameters
constexpr size_t nx = 20; //< Number of grid cells in x-direction
constexpr size_t ny = 20; //< Number of grid cells in y-direction
constexpr double lx = 1.; //< Domain size in x-direction
constexpr double ly = 1.; //< Domain size in y-direction
constexpr double xOrigin = 0.; //< Domain origin in x-direction
constexpr double yOrigin = ly; //< Domain origin in y-direction
constexpr double dx = lx / double(nx); //< Mesh size in x-direction
constexpr double dy = ly / double(ny); //< Mesh size in x-direction

constexpr size_t dimension = 2; //< Space dimension is 2

constexpr double thermalDiff = 1. / 224.36; //< Thermal diffusivity

// Temperature that should be on the interface between the two solvers
// We mimic the Dirichlet boundary from the previous exercise so it should be 1
constexpr double targetTOnInterface = 1.;


static const std::string solverName = "SolidSolverDummy";
static const std::string meshName = "SolidSolverDummyMesh";

// Shorthand definition ofr preCICE constants
// Read cowid  (= co + w + i + d) as
// constant "write initial data"
static const std::string& cowid = precice::constants::actionWriteInitialData();
static const std::string& coric = precice::constants::actionReadIterationCheckpoint();
static const std::string& cowic = precice::constants::actionWriteIterationCheckpoint();


std::array<double, nx*dimension> computeInterfaceCoordinates()
{
  std::array<double, nx*dimension> coordinates;

  for (size_t i = 0; i < nx; ++i) {
    //Set x- and y-coordinate
    const auto x = xOrigin + 0.5 * dx + i * dx;
    const auto y = yOrigin;
    coordinates[i*dimension]   = x;
    coordinates[i*dimension+1] = y;
#ifndef NDEBUG
    std::cout << "Compute coupling point " << i << " at (x,y) = (" << x << ", " << y << ")" << std::endl;
#endif
  }

  return coordinates;
}

// q = - kappa * dT/dy
void computeHeatFluxes( const std::array<double, nx>& temperatureInDomain,
                        const std::array<double, nx>& temperatureInBC,
                        std::array<double, nx>& heatFlux )
{
  // Normal vector of surface is (0,-1)
  const double normalVectorY = -1; //< Normalvector points downwards
  //Compute heat fluxes
  for (size_t i = 0; i < nx; ++i) {
    heatFlux[i] = -thermalDiff * (temperatureInDomain[i] - temperatureInBC[i]) / dy;
    heatFlux[i] *= normalVectorY;
#ifndef NDEBUG
    // Heat flux should be positive (flux out of the domain)
    std::cout << "Heatflux written q[" << i << "] = " << heatFlux[i] << std::endl;
#endif
  }
}

void setBoundaryTemperature( const std::array<double, nx>& temperatureOnInterface,
                             std::array<double, nx>& temperatureInBC )
{
  for (size_t i = 0; i < nx; ++i) {
#ifndef NDEBUG
    std::cout << "Temperature on interface T[" << i << "] = " << temperatureOnInterface[i] << std::endl;
#endif
   //temperatureInBC[i] =  0.5 * (2.*temperatureOnInterface[i] - temperatureInDomain[i]);
    temperatureInBC[i] +=  (temperatureOnInterface[i] - targetTOnInterface);
  }
}

void updateTemperatureInDomain( const std::array<double, nx>& temperatureOnInterface,
                                std::array<double, nx>& temperatureInDomain )
{
  //Compute heat fluxes
  for (size_t i = 0; i < nx; ++i) {
    //temperatureInDomain[i] = 2 - 0.5 * temperatureInBC[i];
    temperatureInDomain[i] -= (temperatureOnInterface[i] - targetTOnInterface );
#ifndef NDEBUG
    std::cout << "Temperature in domain updated to T[" << i << "] = " << temperatureInDomain[i] << std::endl;
#endif
  }
}

// Helper function to create and initialize an array
std::array<double, nx> createAndInitializeArray( const double initialValue = 0. )
{
  std::array<double, nx> arr;
  std::fill( arr.begin(), arr.end(), initialValue );
  return arr;
}

int main( int argc, char *argv[] )
{
  std::string preciceConfigFileName = "precice-config.xml";
  if ( argc == 2 )
  {
    preciceConfigFileName = argv[1];
  }

  // temperatureInDomain is equivalent to the temperature in the
  // boundary/ghost cell of the flow solver
  std::array<double, nx> temperatureInDomain = createAndInitializeArray( 2.0 );
//  std::fill( temperatureInDomain.begin(), temperatureInDomain.end(), 2.0 );

  // temperatureInBC is equivalent to the temperature in the
  // upper most cell of the flow solver that is not a boundary/ghost cell
  std::array<double, nx> temperatureInBC = createAndInitializeArray( 0.0 ) ;
//  std::fill( temperatureInBC.begin(), temperatureInBC.end(), 0.0 );

  //Initialize preCICE
  precice::SolverInterface interface( solverName, 0, 1 );
  interface.configure( preciceConfigFileName );

  // Announce mesh to preCICE
  const int dim = interface.getDimensions();
  assert( dim == dimension );
  const int meshId = interface.getMeshID( meshName );
  std::vector<int> vertexIds(nx, 0);

  // Announce mesh vertices to preCICE
  {
    const std::array<double, nx*dimension> coordinates = computeInterfaceCoordinates();
    interface.setMeshVertices(meshId, int(vertexIds.size()), coordinates.data(), vertexIds.data() );
  }

  // Get data ids
  const int temperatureId = interface.getDataID( "Temperature", meshId );
  const int heatFluxId = interface.getDataID( "Heat-Flux", meshId );

#ifndef NDEBUG
    std::cout << "Data Id temperature: " << temperatureId << std::endl;
    std::cout << "Data Id heat flux: " << heatFluxId << std::endl;
#endif

  std::array<double, nx> temperature = createAndInitializeArray( 0.0 );
  std::array<double, nx> heatFlux = createAndInitializeArray( 0.0 );

  double precice_dt = interface.initialize();

  if ( interface.isActionRequired( cowid ) )
  {
    computeHeatFluxes( temperatureInDomain, temperatureInBC, heatFlux );
    //Write heat flux
    interface.writeBlockScalarData( heatFluxId,
                                    int(vertexIds.size()),
                                    vertexIds.data(),
                                    heatFlux.data() );
    interface.fulfilledAction( cowid );
  }

  interface.initializeData();

  // Starting the time loop (controlled by preCICE)
  while (interface.isCouplingOngoing())
  {
    if ( interface.isActionRequired( cowic ))
    {
      //Save old state
      //I don't do anything. Implicit coupling is not supported here
#ifndef NDEBUG
      std::cout << "HeatDummy: Writing checkpoint!" << std::endl;
#endif
      interface.fulfilledAction( cowic );
    }
    // Read temperature
    interface.readBlockScalarData( temperatureId,
                                   int(vertexIds.size()),
                                   vertexIds.data(),
                                   temperature.data() );

    // Solve problem
    // There is no "real" problem to solve here, but I can update temperatures
    {
      //Set temperature in boundary/ghost cell
      setBoundaryTemperature( temperature, temperatureInBC );
      //Set temperature in interior cell
      updateTemperatureInDomain( temperature, temperatureInDomain );
    }

    //Compute heat fluxes
    computeHeatFluxes( temperatureInDomain, temperatureInBC, heatFlux );

    //Write heat flux
    interface.writeBlockScalarData( heatFluxId,
                                    int(vertexIds.size()),
                                    vertexIds.data(),
                                    heatFlux.data() );
    // Advance time step
    precice_dt = interface.advance( precice_dt );


    if ( interface.isActionRequired( coric ))
    {
      //Load old state
      //I don't do anything. Implicit coupling is not supported here
#ifndef NDEBUG
      std::cout << "HeatDummy: Reading checkpoint!" << std::endl;
#endif
      interface.fulfilledAction( coric );
    }
    else
    {
      //Time step converged
#ifndef NDEBUG
      std::cout << "HeatDummy: Advancing int time!" << std::endl;
#endif
    }
  } // Time loop end

  //Clean up interface/preCICE
  interface.finalize();

#ifndef NDEBUG
      std::cout << "HeatDummy: Terminating!" << std::endl;
#endif

  return 0;
}
