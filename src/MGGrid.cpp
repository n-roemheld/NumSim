#include "MGGrid.h"

// Assumption: nCells are c*2^n_grids

// Use to restrict Staggered Grid to MGGrid on finest level
MGGrid::MGGrid(std::array< int, 2> nCells, std::array< double, 2> meshWidth, FieldVariable p, FieldVariable rhs)
: nCells_(nCells), meshWidth_(meshWidth), p_(p), rhs_(rhs),
  resVec_({nCells[0]+2, nCells[1]+2}, {-.5*meshWidth_[0], -.5*meshWidth_[1]}, meshWidth_)
{};

// Takes values of a finer grid and returns the next-coarser mggrid
MGGrid::MGGrid(std::array< int, 2> nCells, std::array< double, 2> meshWidth)
: nCells_({nCells[0]/2, nCells[1]/2}), meshWidth_({meshWidth[0]*2, meshWidth[1]*2}),
  p_({nCells[0]+2, nCells[1]+2}, {-.5*meshWidth_[0], -.5*meshWidth_[1]}, meshWidth_),
  resVec_({nCells[0]+2, nCells[1]+2}, {-.5*meshWidth_[0], -.5*meshWidth_[1]}, meshWidth_),
  rhs_({nCells[0]+2, nCells[1]+2}, {-.5*meshWidth_[0], -.5*meshWidth_[1]}, meshWidth_)
{};

const std::array< double, 2 > MGGrid::meshWidth() const
{
  return meshWidth_;
};

const std::array< int, 2 > 	MGGrid::nCells() const
{
  return nCells_;
};

const FieldVariable& MGGrid::p() const
{
  return p_;
};

const FieldVariable& MGGrid::rhs() const
{
	return rhs_;
};

const FieldVariable& MGGrid::resVec() const
{
    return resVec_;
};

double MGGrid::p(int i, int j) const
{
  return MGGrid::p_(i,j);
};

double& MGGrid::p(int i, int j)
{
  return MGGrid::p_(i,j);
};

double& MGGrid::rhs(int i, int j)
{
  return MGGrid::rhs_(i,j);
};

double& MGGrid::resVec(int i, int j)
{
    return MGGrid::resVec_(i,j);
};

int MGGrid::pIBegin() const
{
	return 1;
};

int MGGrid::pIEnd() const
{
	return nCells_[0] +1;
};

int MGGrid::pJBegin() const
{
	return 1;
};

int MGGrid::pJEnd() const
{
	return nCells_[1] + 1;
};
