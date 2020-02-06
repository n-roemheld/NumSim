#pragma once
#include "Smoother.h"

class SmootherGaussSeidelUnique : public Smoother
{
public:
    SmootherGaussSeidelUnique(int numberOfIterationsPre, int numberOfIterationsPost);

    void presmooth(std::shared_ptr<MGGrid> mgg);

    void postsmooth(std::shared_ptr<MGGrid> mgg);

protected:
    void smooth(std::shared_ptr<MGGrid> mgg, int numberOfIterations);
};

// GaussSeidel but with setting one node value to 0 to achieve uniqueness of the Neumann poisson problem
