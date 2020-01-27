#include "MGGrid.h"

// Assumption: nCells are c*2^n_grids

// Use to restrict Staggered Grid to MGGrid on finest level
MGGrid::MGGrid(std::array< int, 2> nCells, std::array< double, 2> meshWidth, std::shared_ptr<FieldVariable> p, std::shared_ptr<FieldVariable> rhs)
: nCells_(nCells), meshWidth_(meshWidth), p_(p), rhs_(rhs),
  resVec_({nCells[0]+2, nCells[1]+2}, {-.5*meshWidth_[0], -.5*meshWidth_[1]}, meshWidth_)
{};

// Takes values of a finer grid and returns the next-coarser mggrid
MGGrid::MGGrid(std::array< int, 2> nCells, std::array< double, 2> meshWidth, std::shared_ptr<FieldVariable> rhs)
: nCells_({nCells[0]/2, nCells[1]/2}), meshWidth_({meshWidth[0]*2, meshWidth[1]*2}),
  p_({nCells[0]+2, nCells[1]+2}, {-.5*meshWidth_[0], -.5*meshWidth_[1]}, meshWidth_),
  resVec_({nCells[0]+2, nCells[1]+2}, {-.5*meshWidth_[0], -.5*meshWidth_[1]}, meshWidth_),
  rhs_({nCells[0]+2, nCells[1]+2}, {-.5*meshWidth_[0], -.5*meshWidth_[1]}, meshWidth_)
{};