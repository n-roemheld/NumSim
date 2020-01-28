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
#include "EndSolver.h"
#include "Coarser.h"
#include "CoarserDefault.h"
#include "Multigrid.h"

class Computation
{
public:
	void initialize (int argc, char *argv[]);
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
