#include "Multigrid.h"

Multigrid::Multigrid (std::shared_ptr<Discretization> discretization, double epsilon, int maximumNumberOfIterations, std::unique_ptr<Smoother> sm, std::unique_ptr<Coarser> coa, std::unique_ptr<PressureSolver> es, Cycle cycle) :
    PressureSolver(discretization, epsilon, maximumNumberOfIterations), sm_(sm), coa_(coa), es_(es), cycle_(cycle)
{};


void Multigrid::solve()
{
    std::shared_ptr<MGGrid> mgg = std::make_shared<MGGrid>(discretization_->nCells, discretization_->meshWidth, discretization_->p, discretization_->rhs);
    if(cycle_.recursive)
    {
        MGCycle(cycle_.maxLevel, mgg);
    }
}

Multigrid::computeResVec(std::shared_ptr < MGGrid> mgg)
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

