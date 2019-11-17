#pragma once

#include "PressureSolver.h"
#include "DonorCell.h"
#include "CentralDifferences.h"
#include "SOR.h"
#include "GaussSeidel.h"
#include "SORRedBlack.h"
#include "GaussSeidelRedBlack.h"
#include "output_writer/output_writer_paraview.h"
#include "output_writer/output_writer_text.h"
#include "settings.h"

class Computation
{
public:
	virtual void initialize (int argc, char *argv[]);
	virtual void runSimulation ();

protected:
	Settings settings_;
	std::array< double, 2 > meshWidth_;
	double dt_;

	virtual void computeTimeStepWidth ();
	virtual void applyBoundaryValues ();
	virtual void computePreliminaryVelocities ();
	void computeRightHandSide ();
	void computePressure ();
	virtual void computeVelocities ();

	std::shared_ptr< Discretization > discretization_;
	std::unique_ptr< PressureSolver > pressureSolver_;
	std::unique_ptr< OutputWriterParaview > outputWriterParaview_;
	std::unique_ptr< OutputWriterText > outputWriterText_;


};
