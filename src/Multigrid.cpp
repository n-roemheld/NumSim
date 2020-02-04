#include "Multigrid.h"
#include "GaussSeidel.h"

Multigrid:: Multigrid(std::shared_ptr<Discretization> discretization, std::shared_ptr<Smoother> sm, std::shared_ptr<Coarser> coa, std::shared_ptr<EndSolver> es, Cycle cycle) :
    PressureSolver(discretization, 0, 0), smoother_(sm), coarser_(coa), endSolver_(es), cycle_(cycle) // epsilon und maximumNumberOfIterations einfach auf 0 gesetzt, da nicht gebraucht
{};

void Multigrid::solve()
{
  // std::unique_ptr<GaussSeidel> pressureSolver_ = std::make_unique<GaussSeidel>(discretization_, epsilon_,
  //  maximumNumberOfIterations_);
  //  // discretization_->p(5,5) = -0.7;
  //
  //  pressureSolver_->solve();
   // discretization_->p(5,5) = -7;


    std::shared_ptr<MGGrid> mgg = std::make_shared<MGGrid>(discretization_->nCells(), discretization_->meshWidth(), discretization_->p(), discretization_->rhs());
    if(cycle_.recursive)
    {
      // std::cout << "maxLevel" << cycle_.maxLevel << std::endl;
        MGCycle(cycle_.maxLevel, mgg);
        //neue Werte von mgg in discretization_ schreiben??
    }
    else
    {
        // TODO: iterative
        MGLoop(cycle_.maxLevel, mgg);
    }

    for (int j = discretization_->pJBegin(); j < discretization_->pJEnd(); j++)
    {
      for (int i = discretization_->pIBegin(); i < discretization_->pIEnd(); i++)
      {
        discretization_->p(i,j) = mgg->p().operator()(i,j);
        discretization_->rhs(i,j) = mgg->rhs().operator()(i,j);
      }
    }

};

void Multigrid::MGCycle(int level, std::shared_ptr<MGGrid> mgg)
{
    // std::cout << "Level1: " << level << std::endl;
    if(level == 0)
    {
        endSolver_->solve(mgg);
    }
    else
    {
        smoother_->presmooth(mgg);
        computeResVec(mgg);
        std::shared_ptr<MGGrid> mggc = std::make_shared<MGGrid>(mgg->nCells(), mgg->meshWidth());
        coarser_->restrict(mgg, mggc); // mggCoarse is set complete, also p
        for(int i = 0; i < cycle_.gamma.at(level); i++)
        {
            MGCycle(level-1, mggc);
        }
        coarser_->interpolate(mggc, mgg);
        for(int j = mgg->pJBegin(); j < mgg->pJEnd(); j++)
        {
            for(int i = mgg->pIBegin(); i < mgg->pIEnd(); i++)
            {
                mgg->p(i,j) = mgg->p(i,j) +0.5 * mgg->resVec(i,j); // adding theta ???
            }
        }
        smoother_->postsmooth(mgg);
    }
    //  std::cout << "Level2: " << level << std::endl;
};

void Multigrid::MGLoop(int maxLevel, std::shared_ptr<MGGrid> mgg)
{
    std::vector<std::shared_ptr<MGGrid>> grids;
    grids.push_back(mgg);
    for (int l = 0; l < maxLevel-1; l++)
    {
        smoother_->presmooth(grids.at(l));
        computeResVec(grids.at(l));
        std::shared_ptr<MGGrid> mggc = std::make_shared<MGGrid>(grids.at(l)->nCells(), grids.at(l)->meshWidth());
        coarser_->restrict(mgg, mggc); // mggCoarse is set complete, also p
        grids.push_back(mggc);
    }
    endSolver_->solve(mgg);
    for (int l = maxLevel-2; l >= 0; l--)
    {
        coarser_->interpolate(grids.at(l-1), grids.at(l));
        for(int j = grids.at(l)->pJBegin(); j < grids.at(l)->pJEnd(); j++)
        {
            for(int i = grids.at(l)->pIBegin(); i < grids.at(l)->pIEnd(); i++)
            {
                grids.at(l)->p(i,j) = grids.at(l)->p(i,j) + grids.at(l)->resVec(i,j); // adding theta ???
            }
        }
        smoother_->postsmooth(grids.at(l));
    }
}


void Multigrid::computeResVec(std::shared_ptr < MGGrid> mgg)
{
    std::array<double,2> mW = mgg->meshWidth();
	double dx = mW[0];
	double dy = mW[1];
    for(int j = mgg->pJBegin(); j < mgg->pJEnd(); j++)
    {
        for(int i = mgg->pIBegin(); i < mgg->pIEnd(); i++)
        {
            mgg->resVec(i,j) = mgg->rhs(i,j) - (mgg->p(i+1,j)-2*mgg->p(i,j)+mgg->p(i-1,j))/(dx*dx) - (mgg->p(i,j+1)-2*mgg->p(i,j)+mgg->p(i,j-1))/(dy*dy);
        }
    }
};

void Multigrid::setBoundaryValuesMGGrid(std::shared_ptr<MGGrid> mgg)
{
    // lower p ghost layer
		int j = mgg->pJBegin()-1;
		for(int i = mgg->pIBegin(); i < mgg->pIEnd(); i++)
		{
			mgg->p(i,j) = mgg->p(i,j+1);
		};
		// upper p ghost layer
		j = mgg->pJEnd();
		for(int i = mgg->pIBegin(); i < mgg->pIEnd(); i++)
		{
			mgg->p(i,j) = mgg->p(i,j-1);
		};
		// left p ghost layer
		int i = mgg->pIBegin()-1;
		for(int j = mgg->pJBegin()-1; j < mgg->pJEnd()+1; j++)
		{
			mgg->p(i,j) = mgg->p(i+1,j);
		}
		// right p ghost layer
		i = mgg->pIEnd();
		for(int j = mgg->pJBegin()-1; j < mgg->pJEnd()+1; j++)
		{
			mgg->p(i,j) = mgg->p(i-1,j);
		}
};
