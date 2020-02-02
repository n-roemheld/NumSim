#pragma once

#include "EndSolver.h"

class ESGaussSeidel : public EndSolver
{
  public:
    ESGaussSeidel (double epsilon, int maximumNumberOfIterations);

    void solve(std::shared_ptr<MGGrid> mgg);
};
