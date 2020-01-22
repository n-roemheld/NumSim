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
    Multigrid(std::shared_ptr<Discretization> discretization, double epsilon, int maximumNumberOfIterations, std::unique_ptr<Smoother> sm, std::unique_ptr<Coarser> coa, std::unique_ptr<PressureSolver> es, Cycle cycle);

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
    std::unique_ptr<Smoother> sm_; // smoother
    std::unique_ptr<Coarser> coa_; //coarsening operator
    std::unique_ptr<EndSolver> es_; //endSolver at coarsest level
    Cycle cycle_; // defines cycle structure

};