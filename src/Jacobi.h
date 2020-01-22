#pragma once
#include "Smoother.h"

class Jacobi : public Smoother
{
public: 
    Jacobi(int numberOfIterationsPre, int numberOfIterationsPost);

    // presmooths the mgg with numberOfIterationsPre times Jacobi 
    void presmooth(std::shared_ptr<MGGrid> mgg);

    // postsmooths the mgg with numberOfIterationsPost times Jacobi 
    void postsmooth(std::shared_ptr<MGGrid> mgg);
};