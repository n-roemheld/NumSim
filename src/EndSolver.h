#pragma once

#include "MGGrid.h"
#include <memory>

// does the same as a PressureSolver with a MGGrid instead of a discretization

class EndSolver
{

public:

    EndSolver (double epsilon, int maximumNumberOfIterations);

    virtual void solve(std::shared_ptr<MGGrid> mgg) = 0;

protected:

    void setBoundaryValues(std::shared_ptr<MGGrid> mgg);

    double epsilon_;
    int maximumNumberOfIterations_;

    double compute_res(std::shared_ptr<MGGrid> mgg);
};
