

#include "Computation.h"
#include <math.h>
#include <ctime>

void Computation::initialize (int argc, char *argv[])
{
	settings_.loadFromFile(argv[1]);
	settings_.printSettings();

	//computing meshWidth
	double dx = settings_.physicalSize[0]/settings_.nCells[0];
	double dy = settings_.physicalSize[1]/settings_.nCells[1];
	meshWidth_ = {dx, dy};

	//select DonorCell or CentralDifferences
	if (settings_.useDonorCell == true)
	{
		discretization_ = std::make_shared<DonorCell>(settings_.nCells, meshWidth_, settings_.geometryPVString_, settings_.geometryPV1_, settings_.geometryPV2_, settings_.geometryTString_, settings_.geometryT1_, settings_.alpha, settings_.gamma);
	}
	else
	{
		discretization_ = std::make_shared<CentralDifferences>(settings_.nCells, meshWidth_, settings_.geometryPVString_, settings_.geometryPV1_, settings_.geometryPV2_, settings_.geometryTString_, settings_.geometryT1_);
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

};

void Computation::runSimulation ()
{
	double time = 0;
	double lastOutputTime = 0;
	while(time < settings_.endTime)
	{
		applyObstacleValues();
		applyBoundaryValues();


		if(time == 0) outputWriterParaview_->writeFile(time);

		// std::cout << "time" << time << std::endl;
		// compute dt_ and time
		computeTimeStepWidth();
		if(time+dt_>settings_.endTime) dt_ = settings_.endTime - time;
		// std::cout << "time_step" << dt_ << std::endl;

		// if (time >= nextSnapshotTime)
		// {
		// 	// outputWriterParaview_->writeFile(time);
		// 	// outputWriterText_->writeFile(time);
		// 	// outputWriterText_->writePressureFile();
		// }



		// compute f and g
		computePreliminaryVelocities();
		// outputWriterText_->writeFile(time);

		// compute T
		computeTemperature(); // Reihenfolge?

		//compute rhs
		computeRightHandSide();
		//compute p with SOR or GaussSeidel
		computePressure();
		//compute u and v
		computeVelocities();

		time += dt_;

		if (time - lastOutputTime > settings_.outputFileEveryDt - 1e-4)
		{
			outputWriterParaview_->writeFile(time);
			outputWriterText_->writeFile(time); // todo: disable before submission!
			// outputWriterText_->writePressureFile();
			lastOutputTime = time;
		}
	}
	// output data using VTK if we did not do this in the last time step
	if ( std::fabs( time - lastOutputTime ) > 1e-4 )
	{
		outputWriterParaview_->writeFile(time);
		outputWriterText_->writeFile(time); // todo: disable before submission!
		lastOutputTime = time;
	}

};

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
};

void Computation::applyBoundaryValues ()
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

	// u,f setting
	// lower u ghost layer without corners
	j = discretization_->uJBegin()-1;
	for(int i = discretization_->uIBegin(); i < discretization_->uIEnd()-1; i++)
	{
		discretization_->setBoundaryValues_u_f(lower,i,j);
	};
	// upper u ghost layer without corners
	j = discretization_->uJEnd();
	for(int i = discretization_->uIBegin(); i < discretization_->uIEnd()-1; i++)
	{
		discretization_->setBoundaryValues_u_f(upper,i,j);
	};
	// left u ghost layer with corners
	i = discretization_->uIBegin()-1;
	for(int j = discretization_->uJBegin()-1; j < discretization_->uJEnd()+1; j++)
	{
		discretization_->setBoundaryValues_u_f(left,i,j);
	}
	// right u Nathi-not ghost layer with corners
	i = discretization_->uIEnd()-1;
	for(int j = discretization_->uJBegin()-1; j < discretization_->uJEnd()+1; j++)
	{
		discretization_->setBoundaryValues_u_f(right,i,j);
	}

	// v,g setting
	// lower v ghost layer without corners
	j = discretization_->vJBegin()-1;
	for(int i = discretization_->vIBegin(); i < discretization_->vIEnd(); i++)
	{
		discretization_->setBoundaryValues_v_g(lower,i,j);
	};
	// upper v  Nathi-not ghost layer without corners
	j = discretization_->vJEnd()-1;
	for(int i = discretization_->vIBegin(); i < discretization_->vIEnd(); i++)
	{
		discretization_->setBoundaryValues_v_g(upper,i,j);
	};
	// left v ghost layer with corners
	i = discretization_->vIBegin()-1;
	for(int j = discretization_->vJBegin()-1; j < discretization_->vJEnd(); j++)
	{
		discretization_->setBoundaryValues_v_g(left,i,j);
	}
	// right v ghost layer with corners
	i = discretization_->vIEnd();
	for(int j = discretization_->vJBegin()-1; j < discretization_->vJEnd(); j++)
	{
		discretization_->setBoundaryValues_v_g(right,i,j);
	}
};

void Computation::applyObstacleValues()
{
	// u and f
	for(int j = discretization_->uJBegin(); j < discretization_->uJEnd(); j++)
	{
		for(int i = discretization_->uIBegin(); i < discretization_->uIEnd(); i++)
		{
			discretization_->setObstacleValues_u_f(i,j);
		}
	}
	// v and g
	for(int j = discretization_->vJBegin(); j < discretization_->vJEnd(); j++)
	{
		for(int i = discretization_-> vIBegin(); i < discretization_->vIEnd(); i++)
		{
			discretization_->setObstacleValues_v_g(i,j);
		}
	}
	// T
	for(int j = discretization_->pJBegin(); j < discretization_->pJEnd(); j++)
	{
		for(int i = discretization_-> pIBegin(); i < discretization_->pIEnd(); i++)
		{
			discretization_->setObstacleValues_T(i,j);
		}
	}
}

void Computation::applyObstacleValues2()
{
	// u and f
	for(int j = discretization_->uJBegin(); j < discretization_->uJEnd(); j++)
	{
		int jgeom = j-discretization_->uJBegin()+1;
		for(int i = discretization_->uIBegin(); i < discretization_->uIEnd(); i++)
		{
			int igeom = i-discretization_->uIBegin()+1; // todo: double check!!
			if (discretization_->geometryPVString(igeom, jgeom) == -1)
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
			if (discretization_->geometryPVString(igeom, jgeom) == -1)
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
			if (discretization_->geometryPVString(igeom, jgeom) == -1)
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
									+ (1-settings_.beta*(discretization_->T(i,j)+discretization_->T(i+1,j))/2) * settings_.g[0]);
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
									+ (1-settings_.beta*(discretization_->T(i,j)+discretization_->T(i,j+1))/2) * settings_.g[1]);
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
