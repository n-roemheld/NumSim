#pragma once

#include "PressureSolver.h"

#include <memory>
#include <array>
#include <iostream>
#include "settings.h"
#include "MGGrid.h"
#include "Cycle.h"
#include "Smoother.h"
#include "Coarser.h"
#include "EndSolver.h"
#include "SmootherJacobi.h"
#include "CoarserDefault.h"



class Multigrid : public PressureSolver
{
public:
//    Multigrid(std::shared_ptr<Discretization> discretization, Settings settings_, std::string smoother_name, std::string coarser_name, std::string endSolver_name, Cycle cycle);
    Multigrid(std::shared_ptr<Discretization> discretization, std::shared_ptr<Smoother> sm, std::shared_ptr<Coarser> coa, std::shared_ptr<EndSolver> es, Cycle cycle);

    // wrapper method for multigrid, computes p recursiv or iterative
    void solve();

    // recursive method to perform multigrid
    void MGCycle(int level, std::shared_ptr<MGGrid> mgg);

    // loop method to perform multigrid
    void MGLoop(int maxLevel, std::shared_ptr<MGGrid> mgg);

protected:

    // computes residual vector
    void computeResVec(std::shared_ptr<MGGrid> mgg);

    // set boundary values on MGGrid
    void setBoundaryValuesMGGrid (std::shared_ptr<MGGrid> mgg);

private:

    std::shared_ptr<Smoother> smoother_; // smoother
    std::shared_ptr<Coarser> coarser_; //coarsening operator
    std::shared_ptr<EndSolver> endSolver_; //endSolver at coarsest level
    Cycle cycle_; // defines cycle structure


};
