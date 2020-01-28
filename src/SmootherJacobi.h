#pragma once
#include "Smoother.h"

class SmootherJacobi : public Smoother
{
public: 
    SmootherJacobi(int numberOfIterationsPre, int numberOfIterationsPost);

    // presmooths the mgg with numberOfIterationsPre times Jacobi 
    void presmooth(std::shared_ptr<MGGrid> mgg);

    // postsmooths the mgg with numberOfIterationsPost times Jacobi 
    void postsmooth(std::shared_ptr<MGGrid> mgg);
};
