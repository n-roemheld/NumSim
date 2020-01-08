

#include "Computation.h"
#include <math.h>
#include <ctime>
#include <precice/SolverInterface.hpp>

void Computation::initialize (int argc, char *argv[])
{
	std::cout << "initialize: start" << std::endl;
	settings_.loadFromFile(argv[1]);
	std::cout << "initialize: post load from file" << std::endl;
	settings_.printSettings();

	//computing meshWidth
	double dx = settings_.physicalSize[0]/settings_.nCells[0];
	double dy = settings_.physicalSize[1]/settings_.nCells[1];
	meshWidth_ = {dx, dy};

	// preCICE Adapter
	int rank = 0;
	int size = 1;
	std::cout << "initialize: pre adapter" << std::endl;
	Adapter adapter(settings_.participantName, settings_.preciceConfigFile, rank, size, settings_.vertexSize, settings_.readDataName, settings_.writeDataName);
	std::cout << "initialize: post adapter" << std::endl;

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

	//select SOR or GaussSeidel
	if (settings_.pressureSolver == "SOR")
	{
		pressureSolver_ = std::make_unique<SOR>(discretization_, settings_.epsilon,
		 settings_.maximumNumberOfIterations, settings_.omega);
	}
	else
	{
		pressureSolver_ = std::make_unique<GaussSeidel>(discretization_, settings_.epsilon,
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
	std::vector<double> readData;
	// double *readData = new double[vertexSize];
	std::vector<double> writeData;
	// double *writeData = new double[vertexSize];

	while(time < settings_.endTime)
	// while (discretization_->adapter.isCouplingOngoing()) // implicit
	{
		// if (discretization_->adapter.isActionRequired(discretization_->adapter.cowic))
		// {
		// 	saveOldState();
		// 	discretization_->adapter.fulfilledAction(cowic);
		// }
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

		// compute f and g
		computePreliminaryVelocities();

		// outputWriterText_->writeFile(time);


		outputWriterParaview_->writeFile(time);
		outputWriterText_->writeFile(time);

		// compute T
		computeTemperature(); // Reihenfolge?
		set_writeData(writeData);
		discretization_->adapter.writeData(writeData);

		//compute rhs
		computeRightHandSide();
		//compute p with SOR or GaussSeidel
		computePressure();
		//compute u and v
		computeVelocities();

		discretization_->adapter.advance();

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
	discretization_->adapter.finalize();

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
	discretization_->saveOldStateT();
}

void Computation::reloadOldState()
{
	discretization_->reloadOldStateT();
}


void Computation::set_writeData(std::vector<double> & writeData)
{
	for (int v = 0; v < settings_.vertexSize; v++)
	{
		writeData.at(v) = discretization_->T(settings_.vertex_i.at(v),settings_.vertex_j.at(v));
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

	// compute u_max
	for(int j = 0; j < discretization_->u().size()[1]; j++)
		{
			for (int i = 0; i < discretization_->u().size()[0]; i++)
			{
				if (fabs(discretization_->u(i,j))>u_max) u_max = fabs(discretization_->u(i,j));
			};
		};
	// compute v_max
	for(int j = 0; j < discretization_->v().size()[1]; j++)
		{
			for (int i = 0; i < discretization_->v().size()[0]; i++)
			{
				if (fabs(discretization_->v(i,j))>v_max) v_max = fabs(discretization_->v(i,j));
			};
		};

	double max_dt= settings_.maximumDt;
	// compute mesh dependent time step criterion
	double lim = dx*dx*dy*dy/(dx*dx+dy*dy)*Re/2;
	// check whether pressure diffusion criterion is restricting
	if (lim<max_dt) max_dt= lim;

	lim *= prandtl;
	// check whether temperature diffusion criterion is restricting
	if (lim<max_dt) max_dt= lim;

	// check whether momentum criterions are restricting
	if (u_max>0)
	{
		if(dx/u_max<max_dt) max_dt= dx/u_max;
	}
	if (v_max>0)
	{
		if(dy/v_max<max_dt) max_dt= dy/v_max;
	}
	// multiply with security factor
	dt_ = max_dt*settings_.tau;

	// dt_ = std::min(max_dt*settings_.tau, precice_dt);
};

void Computation::applyBoundaryValues (std::vector<double> & readData)
{
	// setting T boundaries without corners
	discretization_->setBoundaryValues_T(readData, settings_.vertexSize, settings_.vertex_i, settings_.vertex_j);

	// u,f setting
	discretization_->setBoundaryValues_u_f();

	// v,g setting
	discretization_->setBoundaryValues_v_g();
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

// only for not boundary values
void Computation::computePreliminaryVelocities ()
{
	double dt = dt_;
	double Re = settings_.re;
	for(int j = discretization_->uJBegin(); j < discretization_->uJEnd(); j++)
	{
		for (int i = discretization_->uIBegin(); i < discretization_->uIEnd()-1; i++)
		{
			// indices in geometry file (shifted by uIBegin and increased by one at the right (u) and upper(v) boundaries)
			int igeom = i-discretization_-> uIBegin()+1; // todo: double check!!
			int jgeom = j-discretization_-> uJBegin()+1;
			if(discretization_->geometryPVString(igeom,jgeom) == -1 && discretization_->geometryPVString(igeom+1,jgeom) == -1)
			{
				discretization_->f(i,j) = discretization_->u(i,j)
					+ dt*(1/Re*(discretization_->computeD2uDx2(i,j) + discretization_->computeD2uDy2(i,j))
						- discretization_->computeDu2Dx(i,j) - discretization_->computeDuvDy(i,j)
						+ (1 - settings_.beta * (discretization_->T(i,j) + discretization_->T(i+1,j)) / 2) * settings_.g[0]);
						// if (i==3 && j == 3) std::cout << "f" << discretization_->f(i,j)
						// 					<< "Du2Dx2" << discretization_->computeDu2Dx(i,j)
						//           << "Du2Dy2" << discretization_->computeD2uDy2(i,j)
						// 					<< "Du2Dx"  << discretization_->computeDu2Dx(i,j)
						// 					<< "DuvDy"  << discretization_->computeDuvDy(i,j) << '\n';

			}
		};
	};
	for(int j = discretization_->vJBegin(); j < discretization_->vJEnd()-1; j++)
	{
		for (int i = discretization_->vIBegin(); i < discretization_->vIEnd(); i++)
		{
			// indices in geometry file (shifted by uIBegin and increased by one at the right (u) and upper(v) boundaries)
			int igeom = i-discretization_-> vIBegin()+1; // todo: double check!!
			int jgeom = j-discretization_-> vJBegin()+1;
			if(discretization_->geometryPVString(igeom, jgeom) == -1 && discretization_->geometryPVString(igeom, jgeom+1) == -1)
			{
				discretization_->g(i,j) = discretization_->v(i,j)
					+ dt*(1/Re*(discretization_->computeD2vDx2(i,j) + discretization_->computeD2vDy2(i,j))
						- discretization_->computeDv2Dy(i,j) - discretization_->computeDuvDx(i,j)
						+ (1 - settings_.beta*(discretization_->T(i,j) + discretization_->T(i,j+1)) / 2) * settings_.g[1]);
				// std::cout << "g" << discretization_->g(i,j) << "Ts" << discretization_->T(i,j) << " " << discretization_->T(i,j+1) << "g"<< settings_.g[1] << '\n';
			}
		};
	};
};

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
				discretization_->rhs(i,j) = 1/dt*((discretization_->f(i,j)-discretization_->f(i-1,j))/meshWidth_[0]
												+(discretization_->g(i,j)-discretization_->g(i,j-1))/meshWidth_[1]);
			}
		};
	};
};

void Computation::computePressure ()
{
	pressureSolver_->solve();
};

// only for non boundary nodes;
// for boundary nodes the boundary condition would not change anything
//    due to the Neumann boundary condition in p
void Computation::computeVelocities ()
{
	double dt = dt_;
	for(int j = discretization_->uJBegin(); j < discretization_->uJEnd(); j++)
	{
		for (int i = discretization_->uIBegin(); i < discretization_->uIEnd()-1; i++)
		{
			// indices in geometry file (shifted by uIBegin and increased by one at the right (u) and upper(v) boundaries)
			int igeom = i-discretization_-> uIBegin()+1; // todo: double check!!
			int jgeom = j-discretization_-> uJBegin()+1;
			if(discretization_->geometryPVString(igeom,jgeom) == -1 && discretization_->geometryPVString(igeom+1,jgeom) == -1)
			{
				discretization_->u(i,j) = discretization_->f(i,j) - dt*discretization_->computeDpDx(i,j);
				if (i==3 && j==3)
				{
					// std::cout << " u" << discretization_->u(i,j) << " f" << discretization_->f(i,j) << " DpDx" << discretization_->computeDpDx(i,j) << '\n';
				}
			}
		};
	};
	for(int j = discretization_->vJBegin(); j < discretization_->vJEnd()-1; j++)
	{
		for (int i = discretization_->vIBegin(); i < discretization_->vIEnd(); i++)
		{
			// indices in geometry file (shifted by uIBegin and increased by one at the right (u) and upper(v) boundaries)
			int igeom = i-discretization_-> vIBegin()+1; // todo: double check!!
			int jgeom = j-discretization_-> vJBegin()+1;
			if(discretization_->geometryPVString(igeom, jgeom) == -1 && discretization_->geometryPVString(igeom, jgeom+1) == -1)
			{
				discretization_->v(i,j) = discretization_->g(i,j) - dt*discretization_->computeDpDy(i,j);
			}
		};
	};
};

void Computation::computeTemperature()
{
	double dt = dt_;
	FieldVariable T_copy( {settings_.nCells[0]+2, settings_.nCells[1]+2},  {-0.5*meshWidth_[0], -0.5*meshWidth_[1]}, meshWidth_);
	for(int j = discretization_->pJBegin(); j < discretization_->pJEnd(); j++)
	{
		for (int i = discretization_->pIBegin(); i < discretization_->pIEnd(); i++)
		{
			// discretization_->T(i,j) = discretization_->T(i,j)
			T_copy(i,j) = discretization_->T(i,j)
				+ dt*(1/(settings_.re * settings_.prandtl) * ( discretization_->computeD2TDx2(i,j) + discretization_->computeD2TDy2(i,j) )
				- discretization_->computeDuTDx(i,j) - discretization_->computeDvTDy(i,j));
		};
	};

	for(int j = discretization_->pJBegin(); j < discretization_->pJEnd(); j++)
	{
		for (int i = discretization_->pIBegin(); i < discretization_->pIEnd(); i++)
		{
			// indices in geometry file (shifted by uIBegin and increased by one at the right (u) and upper(v) boundaries)
			int igeom = i-discretization_-> pIBegin()+1; // todo: double check!!
			int jgeom = j-discretization_-> pJBegin()+1;
			if(discretization_->geometryPVString(igeom, jgeom) == -1)
			{
				discretization_->T(i,j) = T_copy(i,j);
			}
		};
	};

};
