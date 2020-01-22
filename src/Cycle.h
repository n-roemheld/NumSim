#pragma once 

#include <vector>

struct Cycle
{
    bool recursive = true;
    
    int maxLevel;

    std::vector<int> gamma;
}