#pragma once
#include "Smoother.h"

class SmootherDJacobi : public Smoother
{
public:
    SmootherDJacobi(int numberOfIterationsPre, int numberOfIterationsPost);

    // presmooths the mgg with numberOfIterationsPre times Jacobi
    void presmooth(std::shared_ptr<MGGrid> mgg);

    // postsmooths the mgg with numberOfIterationsPost times Jacobi
    void postsmooth(std::shared_ptr<MGGrid> mgg);

protected:
    void smooth(std::shared_ptr<MGGrid> mgg, int numberOfIterations);
};

// damped Jacobi
