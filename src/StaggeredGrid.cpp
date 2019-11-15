#include "StaggeredGrid.h"

StaggeredGrid::StaggeredGrid(std::array<int,2> nCells, std::array<double,2> meshWidth) :
	nCells_(nCells), meshWidth_(meshWidth),
	u_( {nCells[0]+3, nCells[1]+2},  {-1*meshWidth[0],    -0.5*meshWidth[1]}, meshWidth ),
	v_( {nCells[0]+2, nCells[1]+3},  {-0.5*meshWidth[0], -1*meshWidth[1]}, meshWidth ),
	p_( {nCells[0]+2, nCells[1]+2},  {-0.5*meshWidth[0], -0.5*meshWidth[1]}, meshWidth),
	f_( {nCells[0]+3, nCells[1]+2},  {-1*meshWidth[0],    -0.5*meshWidth[1]}, meshWidth ),
	g_( {nCells[0]+2, nCells[1]+3},  {-0.5*meshWidth[0], -1*meshWidth[1]}, meshWidth ),
	rhs_( {nCells[0]+2, nCells[1]+2},{-0.5*meshWidth[0], -0.5*meshWidth[1]}, meshWidth)
{};

const std::array< double, 2 > StaggeredGrid::meshWidth() const
{
  return meshWidth_;
};

const std::array< int, 2 > 	StaggeredGrid::nCells() const
{
  return nCells_;
};

// removed const after function
const FieldVariable& StaggeredGrid::u() const
{
  return u_;
};

// removed const after function
const FieldVariable& StaggeredGrid::v() const
{
  return v_;
};

// removed const after function
const FieldVariable& StaggeredGrid::p() const
{
  return p_;
};

double StaggeredGrid::u(int i, int j) const
{
  return StaggeredGrid::u_(i,j);
};

double StaggeredGrid::v(int i, int j) const
{
  return StaggeredGrid::v_(i,j);
};

double StaggeredGrid::p(int i, int j) const
{
  return StaggeredGrid::p_(i,j);
};

double& StaggeredGrid::u(int i, int j)
{
  return StaggeredGrid::u_(i,j);
};

double& StaggeredGrid::v(int i, int j)
{
  return StaggeredGrid::v_(i,j);
};

double& StaggeredGrid::p(int i, int j)
{
  return StaggeredGrid::p_(i,j);
};

double& StaggeredGrid::rhs(int i, int j)
{
  return StaggeredGrid::rhs_(i,j);
};

double& StaggeredGrid::f(int i, int j)
{
  return StaggeredGrid::f_(i,j);
};

double& StaggeredGrid::g(int i, int j)
{
  return StaggeredGrid::g_(i,j);
};

double StaggeredGrid::dx() const
{
  return meshWidth_[0];
};

double StaggeredGrid::dy() const
{
  return meshWidth_[1];
};

int StaggeredGrid::uIBegin() const
{
	return 2;
};

int StaggeredGrid::uIEnd() const
{
	return nCells_[0]+2; // reduced by one to account for boundary values // or not? // or not not?
};

int StaggeredGrid::uJBegin() const
{
	return 1;
};

int StaggeredGrid::uJEnd() const
{
	return nCells_[1]+1;
};

int StaggeredGrid::vIBegin() const
{
	return 1;
};

int StaggeredGrid::vIEnd() const
{
	return nCells_[0]+1;
};

int StaggeredGrid::vJBegin() const
{
	return 2;
};

int StaggeredGrid::vJEnd() const
{
	return nCells_[1]+2; // reduced by one to account for boundary values // or not? // or not not?
};

int StaggeredGrid::pIBegin() const
{
	return 1;
};

int StaggeredGrid::pIEnd() const
{
	return nCells_[0] +1;
};

int StaggeredGrid::pJBegin() const
{
	return 1;
};

int StaggeredGrid::pJEnd() const
{
	return nCells_[1] + 1;
};

  void set_partitioning(MPI_rank, ranks_neighbors, is_boundary, nCells)
  {
    partitioning_ = new Partitioning(MPI_rank, ranks_neighbors, is_boundary, nCells);
  };

  const Partitioning& partitioning()
  {
	  return partitioning_;
  };

// old non parallel version !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// #include "StaggeredGrid.h"
//
// StaggeredGrid::StaggeredGrid(std::array<int,2> nCells, std::array<double,2> meshWidth) :
// 	nCells_(nCells), meshWidth_(meshWidth),
// 	u_( {nCells[0]+1, nCells[1]+2},  {0*meshWidth[0],    -0.5*meshWidth[1]}, meshWidth ),
// 	v_( {nCells[0]+2, nCells[1]+1},  {-0.5*meshWidth[0], 0*meshWidth[1]}, meshWidth ),
// 	p_( {nCells[0]+2, nCells[1]+2},  {-0.5*meshWidth[0], -0.5*meshWidth[1]}, meshWidth),
// 	f_( {nCells[0]+1, nCells[1]+2},  {0*meshWidth[0],    -0.5*meshWidth[1]}, meshWidth ),
// 	g_( {nCells[0]+2, nCells[1]+1},  {-0.5*meshWidth[0], 0*meshWidth[1]}, meshWidth ),
// 	rhs_( {nCells[0]+2, nCells[1]+2},{-0.5*meshWidth[0], -0.5*meshWidth[1]}, meshWidth)
// {};
//
// const std::array< double, 2 > StaggeredGrid::meshWidth() const
// {
//   return meshWidth_;
// };
//
// const std::array< int, 2 > 	StaggeredGrid::nCells() const
// {
//   return nCells_;
// };
//
// // removed const after function
// const FieldVariable& StaggeredGrid::u() const
// {
//   return u_;
// };
//
// // removed const after function
// const FieldVariable& StaggeredGrid::v() const
// {
//   return v_;
// };
//
// // removed const after function
// const FieldVariable& StaggeredGrid::p() const
// {
//   return p_;
// };
//
// double StaggeredGrid::u(int i, int j) const
// {
//   return StaggeredGrid::u_(i,j);
// };
//
// double StaggeredGrid::v(int i, int j) const
// {
//   return StaggeredGrid::v_(i,j);
// };
//
// double StaggeredGrid::p(int i, int j) const
// {
//   return StaggeredGrid::p_(i,j);
// };
//
// double& StaggeredGrid::u(int i, int j)
// {
//   return StaggeredGrid::u_(i,j);
// };
//
// double& StaggeredGrid::v(int i, int j)
// {
//   return StaggeredGrid::v_(i,j);
// };
//
// double& StaggeredGrid::p(int i, int j)
// {
//   return StaggeredGrid::p_(i,j);
// };
//
// double& StaggeredGrid::rhs(int i, int j)
// {
//   return StaggeredGrid::rhs_(i,j);
// };
//
// double& StaggeredGrid::f(int i, int j)
// {
//   return StaggeredGrid::f_(i,j);
// };
//
// double& StaggeredGrid::g(int i, int j)
// {
//   return StaggeredGrid::g_(i,j);
// };
//
// double StaggeredGrid::dx() const
// {
//   return meshWidth_[0];
// };
//
// double StaggeredGrid::dy() const
// {
//   return meshWidth_[1];
// };
//
// int StaggeredGrid::uIBegin() const
// {
// 	return 1;
// };
//
// int StaggeredGrid::uIEnd() const
// {
// 	return nCells_[0]+1; // reduced by one to account for boundary values // or not? // or not not?
// };
//
// int StaggeredGrid::uJBegin() const
// {
// 	return 1;
// };
//
// int StaggeredGrid::uJEnd() const
// {
// 	return nCells_[1]+1;
// };
//
// int StaggeredGrid::vIBegin() const
// {
// 	return 1;
// };
//
// int StaggeredGrid::vIEnd() const
// {
// 	return nCells_[0]+1;
// };
//
// int StaggeredGrid::vJBegin() const
// {
// 	return 1;
// };
//
// int StaggeredGrid::vJEnd() const
// {
// 	return nCells_[1]+1; // reduced by one to account for boundary values // or not? // or not not?
// };
//
// int StaggeredGrid::pIBegin() const
// {
// 	return 1;
// };
//
// int StaggeredGrid::pIEnd() const
// {
// 	return nCells_[0] +1;
// };
//
// int StaggeredGrid::pJBegin() const
// {
// 	return 1;
// };
//
// int StaggeredGrid::pJEnd() const
// {
// 	return nCells_[1] + 1;
// };
