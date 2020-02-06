#pragma once
#include "Coarser.h"

class Coarser2 : public Coarser
{
    public:
    Coarser2() : Coarser()
    {};

    // restricts the current MGGrid to the coarser MGGrid and sets also nCells and meshWidth
    void restrict(std::shared_ptr<MGGrid> mggf, std::shared_ptr<MGGrid> mggc);

    // interpolates the coarse MGGrid to the finer MGGrid in p
    void interpolate(std::shared_ptr<MGGrid> mggc, std::shared_ptr<MGGrid> mggf);
};

// coarsening strategy for which coarse grid nodes have same position as fine grid nodes
