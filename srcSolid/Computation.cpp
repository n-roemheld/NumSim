

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

	// Initialize preCICE
	int rank = 0;
	int size = 1;
	// Adapter adapter(settings_.participantName, settings_.preciceConfigFile, rank, size, settings_.vertexSize, settings_.readDataName, settings_.writeDataName);
	// discretization_->adapter.initialize(settings_.meshName, settings_.coords);

	//select DonorCell or CentralDifferences
	if (settings_.useDonorCell == true)
	{
		discretization_ = std::make_shared<DonorCell>(settings_.nCells, meshWidth_, settings_.geometryPVString_, settings_.geometryPVOrientation_, settings_.geometryPV1_, settings_.geometryPV2_, settings_.geometryTString_, settings_.geometryT1_, settings_.alpha, settings_.gamma);
	}
	else
	{
		discretization_ = std::make_shared<CentralDifferences>(settings_.nCells, meshWidth_, settings_.geometryPVString_, settings_.geometryPVOrientation_, settings_.geometryPV1_, settings_.geometryPV2_, settings_.geometryTString_, settings_.geometryT1_);
	}
	discretization_->fillIn(settings_.uInit_, settings_.vInit_, settings_.pInit_, settings_.TInit_);


	//select SOR or GaussSeidel

	temperatureSolver_ = std::make_unique<GaussSeidel>(discretization_, settings_.epsilon,
		 settings_.maximumNumberOfIterations);


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
	// creating interface
	precice::SolverInterface precice(settings_.participantName, 0, 1);
	// basic configuration
	precice.configure(settings_.preciceConfigFile);
	// initialization
    int meshID = precice.getMeshID(settings_.meshName);
    std::vector<int> vertexIDs (settings_.vertexSize,0);
	// set Vertex IDs
	precice.setMeshVertices(meshID,settings_.vertexSize,settings_.coords.data(),vertexIDs.data());
	// initialize precice interface
    double precice_dt = precice.initialize();
	precice.initializeData();
	int writeDataID = precice.getDataID(settings_.writeDataName, meshID);
	int readDataID = precice.getDataID(settings_.readDataName, meshID);
	std::cout << "vs settings: " << settings_.vertexSize << ", vs precice: " << precice.getMeshVertexSize(meshID) << std::endl;
	std::vector<double> readData(settings_.vertexSize,0);
	std::vector<double> writeData(settings_.vertexSize,0);

	while(time < settings_.endTime)
	// while (discretization_->adapter.isCouplingOngoing()) // implicit
	{
		// if (discretization_->adapter.isActionRequired(discretization_->adapter.cowic))
		// {
		// 	saveOldState();
		// 	discretization_->adapter.fulfilledAction(cowic);
		// }

		precice.readBlockScalarData(readDataID, settings_.vertexSize, vertexIDs.data(), readData.data());
		std::cout << "run: post read" << std::endl;
		// applyObstacleValues2();
		applyBoundaryValues(readData);


		// if(time == 0) outputWriterParaview_->writeFile(time);
		if(time == 0) outputWriterText_->writeFile(time);

		// std::cout << "time alpha" << time << '\n';

		// std::cout << "time" << time << std::endl;
		// compute dt_ and time
		computeTimeStepWidth();
		dt_ = std::min(dt_, precice_dt);
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
		// compute T
		computeTemperature(); // Reihenfolge?
		set_writeData(writeData);
		precice.writeBlockScalarData(writeDataID, settings_.vertexSize, vertexIDs.data(), writeData.data());


		precice_dt = precice.advance(dt_);

		// if (discretization_->adapter.isActionRequired(coric))
		// {
		// 	reloadOldState();
		// 	discretization_->adapter.fulfilledAction(coric);
		// }
		// else
		{
			time += dt_;

			if (time - lastOutputTime > settings_.outputFileEveryDt - 1e-4)
			{
				outputWriterParaview_->writeFile(time);
				// outputWriterText_->writeFile(time); // todo: disable before submission!
				// outputWriterText_->writePressureFile();
				lastOutputTime = time;
			}
		}
	}
	precice.finalize();

	// output data using VTK if we did not do this in the last time step
	if ( std::fabs( time - lastOutputTime ) > 1e-4 )
	{
		outputWriterParaview_->writeFile(time);
		// outputWriterText_->writeFile(time); // todo: disable before submission!
		lastOutputTime = time;
	}
};

void Computation::saveOldState()
{
	discretization_->saveOldState();
}

void Computation::reloadOldState()
{
	discretization_->reloadOldState();
}


void Computation::set_writeData(std::vector<double> & writeData)
{
	for(int v = 0; v < settings_.vertexSize; v++)
	{
		int i = settings_.vertex_i.at(v);
		int j = settings_.vertex_j.at(v);

		// std::cout << "here!" << std::endl;

		// neighbour indices
		int in = i;
		int jn = j;
		double h = 0;
		double dx = meshWidth_[0];
		double dy = meshWidth_[1];

		// std::cout << "c" << std::endl;
		// std::cout << "osize" << settings_.orientation_.size() << std::endl;


		int orientation = settings_.orientation_.at(v);

		// std::cout << "here2" << std::endl;

		switch (orientation)
		{
			case 0: std::cout << "no orientation assigned!" << std::endl; break;
			case 1: in = i-1; h = dx; break;
			case 2: in = i+1; h = dx; break;
			case 3: jn = j-1; h = dy; break;
			case 4: jn = j+1; h = dy; break;
			case 5:
			case 6:
			case 7:
			case 8: h = 0; break;
			default: std::cout << "unknown orientation" << std::endl; break;
		}

		// std::cout << "here3!" << std::endl;

		writeData.at(v) = (discretization_->T(i,j)+discretization_->T(in,jn))/2;         // gkkhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhh
	}
}

void Computation::computeTimeStepWidth ()
{
	double dx = meshWidth_[0];
	double dy = meshWidth_[1];
	double Re = settings_.re;
	double prandtl = settings_.prandtl;
	double u_max = 0;
	double v_max = 0;

	double max_dt= settings_.maximumDt;
	// compute mesh dependent time step criterion
	double lim = dx*dx*dy*dy/(dx*dx+dy*dy)*Re/2;
	// check whether pressure diffusion criterion is restricting
	if (lim<max_dt) max_dt= lim;

	lim *= prandtl;
	// check whether temperature diffusion criterion is restricting
	if (lim<max_dt) max_dt= lim;

	// multiply with security factor
	dt_ = max_dt*settings_.tau;

	// dt_ = std::min(max_dt*settings_.tau, precice_dt);
};

void Computation::applyBoundaryValues (std::vector<double> & readData)
{
	// todo: double check indices

	// locations:
	int left = 0;
	int right = 1;
	int lower = 2;
	int upper = 3;

	// setting T boundaries without corners (p grid, all ghost cells)
	// lower T
	int j = discretization_->pJBegin()-1;
	for(int i = discretization_->pIBegin(); i < discretization_->pIEnd(); i++)
	{
		discretization_->setBoundaryValues_T(lower,i,j);
	};
	// upper T
	j = discretization_->pJEnd();
	for(int i = discretization_->pIBegin(); i < discretization_->pIEnd(); i++)
	{
		discretization_->setBoundaryValues_T(upper,i,j);
	};
	// left T
	int i = discretization_->pIBegin()-1;
	for(int j = discretization_->pJBegin()-1; j < discretization_->pJEnd()+1; j++)
	{
		discretization_->setBoundaryValues_T(left,i,j);
	}
	// right T
	i = discretization_->pIEnd();
	for(int j = discretization_->pJBegin()-1; j < discretization_->pJEnd()+1; j++)
	{
		discretization_->setBoundaryValues_T(right,i,j);
	}

	for(int v = 0; v < settings_.vertexSize; v++)
	{
		int i = settings_.vertex_i.at(v);
		int j = settings_.vertex_j.at(v);

		// std::cout << "here!" << std::endl;

		// neighbour indices
		int in = i;
		int jn = j;
		double h = 0;
		double dx = meshWidth_[0];
		double dy = meshWidth_[1];

		// std::cout << "c" << std::endl;
		// std::cout << "osize" << settings_.orientation_.size() << std::endl;


		int orientation = settings_.orientation_.at(v);

		// std::cout << "here2" << std::endl;

		switch (orientation)
		{
			case 0: std::cout << "no orientation assigned!" << std::endl; break;
			case 1: in = i-1; h = dx; break;
			case 2: in = i+1; h = dx; break;
			case 3: jn = j-1; h = dy; break;
			case 4: jn = j+1; h = dy; break;
			case 5:
			case 6:
			case 7:
			case 8: h = 0; break;
			default: std::cout << "unknown orientation" << std::endl; break;
		}

		// std::cout << "here3!" << std::endl;


		discretization_->T(i,j) = h*readData.at(v)*settings_.re*settings_.prandtl        //lkjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjj
		                        + discretization_->T(in,jn);
	}
};


// void Computation::applyBoundaryValues (std::vector<double> & readData)
// {
// 	// setting T boundaries without corners
// 	discretization_->setBoundaryValues_T(readData, settings_.vertexSize, settings_.vertex_i, settings_.vertex_j);
//
// 	// u,f setting
// 	discretization_->setBoundaryValues_u_f();
//
// 	// v,g setting
// 	discretization_->setBoundaryValues_v_g();
// };

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
	// u and f
	for(int j = discretization_->uJBegin(); j < discretization_->uJEnd(); j++)
	{
		int jgeom = j-discretization_->uJBegin()+1;
		for(int i = discretization_->uIBegin(); i < discretization_->uIEnd(); i++)
		{
			int igeom = i-discretization_->uIBegin()+1; // todo: double check!!
			if (discretization_->geometryPVString(igeom, jgeom) == 5)
			{
				discretization_->setObstacleValues_u_f2(i,j);
			}
		}
	}
	// v and g
	for(int j = discretization_->vJBegin(); j < discretization_->vJEnd(); j++)
	{
		int jgeom = j-discretization_->vJBegin()+1;
		for(int i = discretization_-> vIBegin(); i < discretization_->vIEnd(); i++)
		{
			int igeom = i-discretization_->vIBegin()+1; // todo: double check!!
			if (discretization_->geometryPVString(igeom, jgeom) == 5)
			{
				discretization_->setObstacleValues_v_g2(i,j);
			}
		}
	}
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
