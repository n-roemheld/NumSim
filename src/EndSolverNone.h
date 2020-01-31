#include "EndSolver.h"

class EndSolverNone : public EndSolver
{
    public:

    EndSolverNone(double epsilon, int maximumNumberOfIterations)
    : EndSolver(epsilon, maximumNumberOfIterations)
    {};

    void solve(std::shared_ptr<MGGrid> mgg)
    {};
};