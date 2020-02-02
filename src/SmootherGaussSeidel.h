#pragma once 
#include "Smoother.h"

class SmootherGaussSeidel : public Smoother
{
public:
    SmootherGaussSeidel(int numberOfIterationsPre, int numberOfIterationsPost);

    void presmooth(std::shared_ptr<MGGrid> mgg);

    void postsmooth(std::shared_ptr<MGGrid> mgg);

protected:
    void smooth(std::shared_ptr<MGGrid> mgg, int numberOfIterations);
};
