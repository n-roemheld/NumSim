#pragma once

#include <memory>
#include "Smoother.h"
#include "Coarser.h"
#include "MGGrid.h"
#include "PressureSolver.h"
#include "Discretization.h"
#include "Cycle.h"
#include "EndSolver.h"



class Multigrid : public PressureSolver
{
public:
    Multigrid(std::shared_ptr<Discretization> discretization, double epsilon, int maximumNumberOfIterations, std::shared_ptr<Smoother> sm, std::shared_ptr<Coarser> coa, std::shared_ptr<EndSolver> es, Cycle cycle);

    // wrapper method for multigrid, computes p recursiv or iterative
    void solve();

    // recursive method to perform multigrid
    void MGCycle(int level, std::shared_ptr<MGGrid> mgg); 

protected:
    
    // computes residual vector
    void computeResVec(std::shared_ptr<MGGrid> mgg);

    // set boundary values on MGGrid
    void setBoundaryValuesMGGrid (std::shared_ptr<MGGrid> mgg);

private:
    std::shared_ptr<Smoother> sm_; // smoother
    std::shared_ptr<Coarser> coa_; //coarsening operator
    std::shared_ptr<EndSolver> es_; //endSolver at coarsest level
    Cycle cycle_; // defines cycle structure

};