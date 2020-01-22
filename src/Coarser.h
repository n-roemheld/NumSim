#pragma once

#include <memory>
#include <vector>
#include "MGGrid.h"

class Coarser
{
public:

    Coarser();

    // restricts the current MGGrid to the coarser MGGrid and sets also nCells and meshWidth 
    virtual std::shared_ptr<MGGrid> restrict(std::shared_ptr<MGGrid> mgg) = 0;

    // interpolates the coarse MGGrid to the finer MGGrid in p
    virtual std::shared_ptr<FieldVariable> interpolate(std::shared_ptr<MGGrid> mggCoarse) = 0;
};