#pragma once

#include "EndSolver.h"

class MGGaussSeidel : public EndSolver
{
  public:
    MGGaussSeidel (double epsilon, int maximumNumberOfIterations);

    void solve(std::shared_ptr<MGGrid> mgg);
}
