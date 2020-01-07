#pragma once

#include "TemperatureSolver.h"
#include "DonorCell.h"
#include "CentralDifferences.h"
#include "SOR.h"
#include "GaussSeidel.h"
#include "output_writer/output_writer_paraview.h"
#include "output_writer/output_writer_text.h"
#include "settings.h"

class Computation
{
public:
	void initialize (int argc, char *argv[]);
	void runSimulation ();
private:
	void computeTimeStepWidth ();
	void applyBoundaryValues (double * readData);
	void computeRightHandSide ();
	void computeTemperature();
	// void applyObstacleValues();
	void applyObstacleValues2();
	void set_writeData(double* writeData);

	Settings settings_;
	std::shared_ptr< Discretization > discretization_;
	std::unique_ptr< TemperatureSolver > temperatureSolver_;
	std::unique_ptr< OutputWriterParaview > outputWriterParaview_;
	std::unique_ptr< OutputWriterText > outputWriterText_;
	std::array< double, 2 > meshWidth_;
	double dt_;

};
