#pragma once

#include <memory>
#include <vector>
#include "MGGrid.h"

class Coarser
{
public:

    Coarser();

    virtual std::vector<double> restrict(std::vector<double> fineResVec, std::shared_ptr<MGGrid> mgg) = 0;

    virtual std::vector<double> interpolate(std::vector<double> coaDefVec, std::shared_ptr<MGGrid> mgg) = 0;
};