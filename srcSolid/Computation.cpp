

#include "Computation.h"
#include <math.h>
#include <ctime>
#include <precice/SolverInterface.hpp>

void Computation::initialize (int argc, char *argv[])
{
	settings_.loadFromFile(argv[1]);
	settings_.printSettings();

	//computing meshWidth
	double dx = settings_.physicalSize[0]/settings_.nCells[0];
	double dy = settings_.physicalSize[1]/settings_.nCells[1];
	meshWidth_ = {dx, dy};

	// preCICE Adapter
	int rank = 0;
	int size = 1;
	Adapter adapter(settings_.participantName, settings_.preciceConfigFile, rank, size, settings_.vertexSize, settings_.readDataName, settings_.writeDataName);

	//select DonorCell or CentralDifferences
	if (settings_.useDonorCell == true)
	{
		discretization_ = std::make_shared<DonorCell>(settings_.nCells, meshWidth_, settings_.geometryPVString_, settings_.geometryPVOrientation_, settings_.geometryPV1_, settings_.geometryPV2_, settings_.geometryTString_, settings_.geometryT1_, settings_.alpha, settings_.gamma, adapter);
	}
	else
	{
		discretization_ = std::make_shared<CentralDifferences>(settings_.nCells, meshWidth_, settings_.geometryPVString_, settings_.geometryPVOrientation_, settings_.geometryPV1_, settings_.geometryPV2_, settings_.geometryTString_, settings_.geometryT1_, adapter);
	}
	discretization_->fillIn(settings_.uInit_, settings_.vInit_, settings_.pInit_, settings_.TInit_);

// always select GaussSeidel
//	//select SOR or GaussSeidel
//	if (settings_.temperatureSolver == "SOR")
//	{
//		temperatureSolver_ = std::make_unique<SOR>(discretization_, settings_.epsilon,
//		 settings_.maximumNumberOfIterations, settings_.omega);
//	}
//	else
	{
		temperatureSolver_ = std::make_unique<GaussSeidel>(discretization_, settings_.epsilon,
		 settings_.maximumNumberOfIterations);
	}

	//initialize outputWriters
 	outputWriterParaview_ = std::make_unique<OutputWriterParaview>(discretization_);
	outputWriterText_ = std::make_unique<OutputWriterText>(discretization_);

	// outputWriterText_->writeFile(-1);
};

void Computation::runSimulation ()
{
	double time = 0;
	double lastOutputTime = 0;

	// Initialize preCICE

	discretization_->adapter.initialize(settings_.meshName, settings_.coords);
	int vertexSize = discretization_->adapter.getVertexSize();
	double *readData = new double[vertexSize];
	double *writeData = new double[vertexSize];

	while(time < settings_.endTime)
	{
		discretization_->adapter.readData(readData);
		applyObstacleValues2();
		applyBoundaryValues(readData);


		// if(time == 0) outputWriterParaview_->writeFile(time);
		if(time == 0) outputWriterText_->writeFile(time);

		// std::cout << "time alpha" << time << '\n';

		// std::cout << "time" << time << std::endl;
		// compute dt_ and time
		computeTimeStepWidth();
		dt_ = discretization_->adapter.get_dt(dt_);
		if(time+dt_>settings_.endTime) dt_ = settings_.endTime - time;
		// std::cout << "time_step" << dt_ << std::endl;

		// if (time - lastOutputTime > settings_.outputFileEveryDt - 1e-4)
		// {
			// outputWriterParaview_->writeFile(time);
			// outputWriterText_->writeFile(time);
			// outputWriterText_->writePressureFile();
		// }


		std::cout << "time1 " << time << '\n';
		std::cout << "dt" << dt_ << '\n';

		// outputWriterText_->writeFile(time);


		outputWriterParaview_->writeFile(time);
		outputWriterText_->writeFile(time);

		//compute rhs
		computeRightHandSide();

		// compute T with GaussSeidel
		computeTemperature(); // Reihenfolge?
		set_writeData(writeData);
		discretization_->adapter.writeData(writeData);

		discretization_->adapter.advance();

		time += dt_;

		if (time - lastOutputTime > settings_.outputFileEveryDt - 1e-4)
		{
			outputWriterParaview_->writeFile(time);
			// outputWriterText_->writeFile(time); // todo: disable before submission!
			// outputWriterText_->writePressureFile();
			lastOutputTime = time;
		}
	}
	discretization_->adapter.finalize();

	// output data using VTK if we did not do this in the last time step
	if ( std::fabs( time - lastOutputTime ) > 1e-4 )
	{
		outputWriterParaview_->writeFile(time);
		// outputWriterText_->writeFile(time); // todo: disable before submission!
		lastOutputTime = time;
	}

};

void Computation::set_writeData(double* writeData)
{
	for (int v = 0; v < settings_.vertexSize; v++)
	{
		writeData[v] = discretization_->T(settings_.vertex_i.at(v),settings_.vertex_j.at(v));
	}
}

void Computation::computeTimeStepWidth ()
{
	double dx = meshWidth_[0];
	double dy = meshWidth_[1];
	double Re = settings_.re;
	double prandtl = settings_.prandtl;

	double max_dt= settings_.maximumDt;
	// compute mesh dependent time step criterion
	double lim = dx*dx*dy*dy/(dx*dx+dy*dy)*Re/2;


	lim *= prandtl;
	// check whether temperature diffusion criterion is restricting
	if (lim<max_dt) max_dt= lim;

	// multiply with security factor
	dt_ = max_dt*settings_.tau;

	// dt_ = std::min(max_dt*settings_.tau, precice_dt);
};

void Computation::applyBoundaryValues (double * readData)
{
	// setting T boundaries without corners
	discretization_->setBoundaryValues_T(readData, settings_.vertexSize, settings_.vertex_i, settings_.vertex_j);
};

// void Computation::applyObstacleValues()
// {
// 	// u and f
// 	for(int j = discretization_->uJBegin(); j < discretization_->uJEnd(); j++)
// 	{
// 		for(int i = discretization_->uIBegin(); i < discretization_->uIEnd(); i++)
// 		{
// 			discretization_->setObstacleValues_u_f(i,j);
// 		}
// 	}
// 	// v and g
// 	for(int j = discretization_->vJBegin(); j < discretization_->vJEnd(); j++)
// 	{
// 		for(int i = discretization_-> vIBegin(); i < discretization_->vIEnd(); i++)
// 		{
// 			discretization_->setObstacleValues_v_g(i,j);
// 		}
// 	}
// 	// T
// 	for(int j = discretization_->pJBegin(); j < discretization_->pJEnd(); j++)
// 	{
// 		for(int i = discretization_-> pIBegin(); i < discretization_->pIEnd(); i++)
// 		{
// 			discretization_->setObstacleValues_T(i,j);
// 		}
// 	}
// }

void Computation::applyObstacleValues2()
{
	// T
	for(int j = discretization_->pJBegin(); j < discretization_->pJEnd(); j++)
	{
		int jgeom = j-discretization_->pJBegin()+1;
		for(int i = discretization_-> pIBegin(); i < discretization_->pIEnd(); i++)
		{
			int igeom = i-discretization_->pIBegin()+1; // todo: double check!!
			if (discretization_->geometryPVString(igeom, jgeom) == 5)
			{
				discretization_->setObstacleValues_T2(i,j);
			}
		}
	}
}

void Computation::computeRightHandSide ()
{
	double dt = dt_;
	for(int j = discretization_->pJBegin(); j < discretization_->pJEnd(); j++)
	{
		for (int i = discretization_->pIBegin(); i < discretization_->pIEnd(); i++)
		{
			// indices in geometry file (shifted by uIBegin and increased by one at the right (u) and upper(v) boundaries)
			int igeom = i-discretization_-> pIBegin()+1; // todo: double check!!
			int jgeom = j-discretization_-> pJBegin()+1;
			if(discretization_->geometryPVString(igeom, jgeom) == -1)
			{
				discretization_->rhs(i,j) = discretization_->T(i,j);
			}
		};
	};
};

void Computation::computeTemperature ()
{
	temperatureSolver_->solve(dt_, settings_.heatDiffusivity);
};

