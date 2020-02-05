#pragma once

#include "PressureSolver.h"
#include "DonorCell.h"
#include "CentralDifferences.h"
#include "SOR.h"
#include "GaussSeidel.h"
#include "output_writer/output_writer_paraview.h"
#include "output_writer/output_writer_text.h"
#include "settings.h"

#include "Smoother.h"
#include "SmootherJacobi.h"
#include "SmootherDJacobi.h"
#include "SmootherGaussSeidel.h"
#include "SmootherGaussSeidelUnique.h"

#include "Multigrid.h"
#include "EndSolver.h"
#include "EndSolverNone.h"
#include "ESGaussSeidel.h"
#include "Coarser.h"
#include "CoarserDefault.h"
#include "CoarserLinear.h"
#include "Coarser2.h"

class Computation
{
public:
	//----changed for creating data--------
	// void initialize (int argc, char *argv[]);
	void initialize (int argc, char *argv[], std::string smootherString, std::string coarserString, std::string endSolverString, int mL, int noIPre, int noIPost, std::string cycleString);
	//-------------------------------------
	void runSimulation ();

private:
	void computeTimeStepWidth ();
	void applyBoundaryValues ();
	void computePreliminaryVelocities ();
	void computeRightHandSide ();
	void computePressure ();
	void computeVelocities ();

	Settings settings_;
	std::shared_ptr< Discretization > discretization_;
	std::unique_ptr< PressureSolver > pressureSolver_;
	std::unique_ptr< OutputWriterParaview > outputWriterParaview_;
	std::unique_ptr< OutputWriterText > outputWriterText_;
	std::array< double, 2 > meshWidth_;
	double dt_;


};
