#pragma once

#include <memory>
#include <vector>
#include "MGGrid.h"

class Coarser
{
public:

    Coarser()
    {};

    // restricts the current MGGrid to the coarser MGGrid and sets also nCells and meshWidth 
    virtual void restrict(std::shared_ptr<MGGrid> mggf, std::shared_ptr<MGGrid> mggc) = 0;

    // interpolates the coarse MGGrid to the finer MGGrid in p
    virtual void interpolate(std::shared_ptr<MGGrid> mggc, std::shared_ptr<MGGrid> mggf) = 0;
};