#include "Multigrid.h"
#include <array>

Multigrid::Multigrid (std::shared_ptr<Discretization> discretization, double epsilon, int maximumNumberOfIterations, std::shared_ptr<Smoother> sm, std::shared_ptr<Coarser> coa, std::shared_ptr<EndSolver> es, Cycle cycle) :
    PressureSolver(discretization, epsilon, maximumNumberOfIterations), sm_(sm), coa_(coa), es_(es), cycle_(cycle)
{};


void Multigrid::solve()
{
    std::shared_ptr<FieldVariable> p = std::make_shared<FieldVariable>(discretization_->p()); // wird eine Kopie erstellt??
    std::shared_ptr<FieldVariable> rhs = std::make_shared<FieldVariable>(discretization_->rhs()); // wird eine Kopie erstellt??
    std::shared_ptr<MGGrid> mgg = std::make_shared<MGGrid>(discretization_->nCells(), discretization_->meshWidth(), p, rhs);
    if(cycle_.recursive)
    {
        MGCycle(cycle_.maxLevel, mgg);
        //neue Werte von mgg in discretization_ schreiben??
    }
    else
    {
        // TODO: iterative
    }

};

void Multigrid::MGCycle(int level, std::shared_ptr<MGGrid> mgg)
{
    if(level == 0)
    {
        es_->solve(mgg);
    }
    else 
    {
        sm_->presmooth(mgg);
        computeResVec(mgg);
        MGGrid mggc_obj = MGGrid(mgg->nCells(), mgg->meshWidth());
        std::shared_ptr<MGGrid> mggc = std::make_shared<MGGrid>(mggc_obj);
        coa_->restrict(mgg, mggc); // mggCoarse is set complete, also p
        for(int i = 0; i < cycle_.gamma[level]; i++)
        {
            MGCycle(level-1, mggc);
        }
        coa_->interpolate(mggc, mgg);
        for(int j = mgg->pJBegin(); j < mgg->pJEnd(); j++)
        {
            for(int i = mgg->pIBegin(); i < mgg->pIEnd(); i++)
            {
                mgg->p(i,j) = mgg->p(i,j) + mgg->resVec(i,j); // adding theta ???
            }
        }
        sm_->postsmooth(mgg);
    }
};



void Multigrid::computeResVec(std::shared_ptr < MGGrid> mgg)
{
    FieldVariable p = mgg->p();
    std::array<double,2> mW = mgg->meshWidth();
	double dx = mW[0];
	double dy = mW[1];
    for(int j = mgg->pJBegin(); j < mgg->pJEnd(); j++)
    {
        for(int i = mgg->pIBegin(); i < mgg->pIEnd(); i++)
        {
            mgg->resVec(i,j) = mgg->rhs(i,j) - (p(i+1,j)-2*p(i,j)+p(i-1,j))/(dx*dx) + (p(i,j+1)-2*p(i,j)+p(i,j-1))/(dy*dy);
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

