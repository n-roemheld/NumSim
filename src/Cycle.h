#pragma once 

#include <vector>

struct Cycle
{
    bool recursive = true; // defines wheather multigrid is done recursive or iterative
    
    int maxLevel; // maxLevel for multigrid

    std::vector<int> gamma; // v.at(level) defines times to perform multigrid on special level
}