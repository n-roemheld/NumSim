#pragma once
#include "Smoother.h"

class Jacobi : public Smoother
{
public: 
    Jacobi(int numberOfIterations);

    void smooth(MGGrid);
};