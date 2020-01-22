#pragma once

#include <memory>
#include "MGGrid.h"

class Smoother
{
public:

    Smoother (int numberOfIterationsPre, int numberOfIterationsPost);

    // performs whole presmoothing depending on numberOfIterationsPre_
    virtual void presmooth(MGGrid mgg) = 0;

    // performs whole postsmoothing depending on numberOfIterationsPost_
    virtual void postsmooth(MGGrid mgg) = 0;

protected:

    int numberOfIterationsPre_;

    int numberOfIterationsPost_;
};