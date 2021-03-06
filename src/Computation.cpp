
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
		discretization_ = std::make_shared<DonorCell>(settings_.nCells, meshWidth_, settings_.alpha);
	}
	else
	{
		discretization_ = std::make_shared<CentralDifferences>(settings_.nCells, meshWidth_);
	}

	//select SOR or GaussSeidel
	if (settings_.pressureSolver == "SOR")
	{
		pressureSolver_ = std::make_unique<SORRedBlack>(discretization_, settings_.epsilon,
		 settings_.maximumNumberOfIterations, settings_.omega);
	}
	else
	{
		pressureSolver_ = std::make_unique<GaussSeidel>(discretization_, settings_.epsilon,
		 settings_.maximumNumberOfIterations);
	}

	//initialize outputWriters
 	// outputWriterParaview_ = std::make_unique<OutputWriterParaview>(discretization_);
	// outputWriterText_ = std::make_unique<OutputWriterText>(discretization_);

};

void Computation::runSimulation ()
{
	double time = 0;
	while(time<settings_.endTime)
	{
		applyBoundaryValues();

		// std::cout << "time" << time << std::endl;
		// compute dt_ and time
		computeTimeStepWidth();
		if(time+dt_>settings_.endTime) dt_ = settings_.endTime - time;
		time += dt_;
		//std::cout << "time_step" << dt_ << std::endl;

		// compute f and g
		computePreliminaryVelocities();
		// outputWriterText_->writeFile(time);

		//compute rhs
		computeRightHandSide();
		//compute p with SOR or GaussSeidel
		computePressure();
		//compute u and v
		computeVelocities();

		outputWriterParaview_->writeFile(time);
		outputWriterText_->writeFile(time);
		outputWriterText_->writePressureFile();
	}
};

void Computation::computeTimeStepWidth ()
{
	double dx = meshWidth_[0];
	double dy = meshWidth_[1];
	double Re = settings_.re;
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
	// check whether mesh dependent time step criterion is restricting
	if (lim<max_dt) max_dt= lim;

	// check whether velocity dependent time step criterions are restricting
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
	// u,f setting
	// lower u ghost layer without corners
	int j = discretization_->uJBegin()-1;
	for(int i = discretization_->uIBegin(); i < discretization_->uIEnd()-1; i++)
	{
		discretization_->u(i,j) = 2*settings_.dirichletBcBottom[0]-discretization_->u(i,j+1);
		discretization_->f(i,j) = discretization_->u(i,j);
	};
	// upper u ghost layer without corners
	j = discretization_->uJEnd();
	for(int i = discretization_->uIBegin(); i < discretization_->uIEnd()-1; i++)
	{
		discretization_->u(i,j) = 2*settings_.dirichletBcTop[0]-discretization_->u(i,j-1);
		discretization_->f(i,j) = discretization_->u(i,j);
	};
	// left u ghost layer with corners
	int i = discretization_->uIBegin()-1;
	for(int j = discretization_->uJBegin()-1; j < discretization_->uJEnd()+1; j++)
	{
		discretization_->u(i,j) = settings_.dirichletBcLeft[0];
		discretization_->f(i,j) = discretization_->u(i,j);
	}
	// right u Nathi-not ghost layer with corners
	i = discretization_->uIEnd()-1;
	for(int j = discretization_->uJBegin()-1; j < discretization_->uJEnd()+1; j++)
	{
		discretization_->u(i,j) = settings_.dirichletBcRight[0];
		discretization_->f(i,j) = discretization_->u(i,j);
	}

	// v,g setting
	// lower v ghost layer without corners
	j = discretization_->vJBegin()-1;
	for(int i = discretization_->vIBegin(); i < discretization_->vIEnd(); i++)
	{
		discretization_->v(i,j) = settings_.dirichletBcBottom[1];
		discretization_->g(i,j) = discretization_->v(i,j);
	};
	// upper v  Nathi-not ghost layer without corners
	j = discretization_->vJEnd()-1;
	for(int i = discretization_->vIBegin(); i < discretization_->vIEnd(); i++)
	{
		discretization_->v(i,j) = settings_.dirichletBcTop[1];
		discretization_->g(i,j) = discretization_->v(i,j);
	};
	// left v ghost layer with corners
	i = discretization_->vIBegin()-1;
	for(int j = discretization_->vJBegin()-1; j < discretization_->vJEnd(); j++)
	{
		discretization_->v(i,j) = 2*settings_.dirichletBcLeft[1]-discretization_->v(i+1,j);
		discretization_->g(i,j) = discretization_->v(i,j);
	}
	// right v ghost layer with corners
	i = discretization_->vIEnd();
	for(int j = discretization_->vJBegin()-1; j < discretization_->vJEnd(); j++)
	{
		discretization_->v(i,j) = 2*settings_.dirichletBcRight[1]-discretization_->v(i-1,j);
		discretization_->g(i,j) = discretization_->v(i,j);
	}
};

// only for not boundary values
void Computation::computePreliminaryVelocities ()
{
	double dt = dt_;
	double Re = settings_.re;
	for(int j = discretization_->uJBegin(); j < discretization_->uJEnd(); j++)
				{
					for (int i = discretization_->uIBegin(); i < discretization_->uIEnd()-1; i++)
					{
						discretization_->f(i,j) = discretization_->u(i,j)
								+ dt*(1/Re*(discretization_->computeD2uDx2(i,j) + discretization_->computeD2uDy2(i,j))
									- discretization_->computeDu2Dx(i,j) - discretization_->computeDuvDy(i,j) + settings_.g[0]);
					};
				};
	for(int j = discretization_->vJBegin(); j < discretization_->vJEnd()-1; j++)
				{
					for (int i = discretization_->vIBegin(); i < discretization_->vIEnd(); i++)
					{
						discretization_->g(i,j) = discretization_->v(i,j)
								+ dt*(1/Re*(discretization_->computeD2vDx2(i,j) + discretization_->computeD2vDy2(i,j))
									- discretization_->computeDv2Dy(i,j) - discretization_->computeDuvDx(i,j) + settings_.g[1]);
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
					// old version
					// discretization_->rhs(i,j) = 1/dt*((discretization_->f(i,j)-discretization_->f(i-1,j))/meshWidth_[0]
					// 							     +(discretization_->g(i,j)-discretization_->g(i,j-1))/meshWidth_[1]);
	  			 // new version
					 discretization_->rhs(i,j) = 1/dt*((discretization_->f(i+1,j)-discretization_->f(i,j))/meshWidth_[0]
					 												     +(discretization_->g(i,j+1)-discretization_->g(i,j))/meshWidth_[1]);

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
					// old version
					// discretization_->u(i,j) = discretization_->f(i,j) - dt*discretization_->computeDpDx(i,j);
					// new version
					discretization_->u(i,j) = discretization_->f(i,j) - dt*discretization_->computeDpDx(i-1,j);
				};
			};
	for(int j = discretization_->vJBegin(); j < discretization_->vJEnd()-1; j++)
			{
				for (int i = discretization_->vIBegin(); i < discretization_->vIEnd(); i++)
				{
					// old version
					// discretization_->v(i,j) = discretization_->g(i,j) - dt*discretization_->computeDpDy(i,j);
					// new version
					discretization_->v(i,j) = discretization_->g(i,j) - dt*discretization_->computeDpDy(i,j-1);
				};
			};
};
