#pragma once

#include <memory>
#include "MGGrid.h"

class Smoother
{
public:

    Smoother (int numberOfIterationsPre, int numberOfIterationsPost);

    // performs whole presmoothing depending on numberOfIterationsPre_
    virtual void presmooth(std::shared_ptr<MGGrid> mgg) = 0;

    // performs whole postsmoothing depending on numberOfIterationsPost_
    virtual void postsmooth(std::shared_ptr<MGGrid> mgg) = 0;

protected:

    void setBoundaryValues(std::shared_ptr<MGGrid> mgg);

    int numberOfIterationsPre_; // number of Iterations for presmoothing

    int numberOfIterationsPost_; // number of iterations for postsmoothing
};