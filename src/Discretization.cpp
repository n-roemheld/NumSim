#include "Discretization.h"

Discretization::Discretization(std::array<int, 2> nCells,
		std::array<double, 2> meshWidth) :
		StaggeredGrid(nCells, meshWidth)
{};

double Discretization::computeD2uDx2(int i, int j) const
{
	return (u_(i + 1, j) - 2 * u_(i, j) + u_(i - 1, j))
			/ (meshWidth_[0] * meshWidth_[0]);
};

double Discretization::computeD2uDy2(int i, int j) const
{
	return (u_(i, j + 1) - 2 * u_(i, j) + u_(i, j - 1))
			/ (meshWidth_[1] * meshWidth_[1]);
};

double Discretization::computeD2vDx2(int i, int j) const
{
	return (v_(i + 1, j) - 2 * v_(i, j) + v_(i - 1, j))
			/ (meshWidth_[0] * meshWidth_[0]);
};

double Discretization::computeD2vDy2(int i, int j) const
{
	return (v_(i, j + 1) - 2 * v_(i, j) + v_(i, j - 1))
			/ (meshWidth_[1] * meshWidth_[1]);
};

double Discretization::computeDpDx(int i, int j) const
{
	return (p_(i+1,j)-p_(i,j))/meshWidth_[0];
};

double Discretization::computeDpDy(int i, int j) const
{
	return (p_(i,j+1)-p_(i,j))/meshWidth_[1];
};

double Discretization::computeD2TDx2(int i, int j) const
{
	return (T_(i+1,j)-2*T_(i,j)+T_(i-1,j))/meshWidth_[0];
};

double Discretization::computeD2TDy2(int i, int j) const
{
	return (T_(i,j+1)-2*T_(i,j)+T_(i,j-1))/meshWidth_[1];
};
