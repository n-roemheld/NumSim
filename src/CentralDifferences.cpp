#include "CentralDifferences.h"


CentralDifferences::CentralDifferences(std::array< int, 2 > nCells, std::array< double, 2 > meshWidth) :
	Discretization(nCells, meshWidth)
{};

CentralDifferences::CentralDifferences(std::array< int, 2 > nCells, std::array< double, 2 > meshWidth, Partitioning parti) :
	Discretization(nCells, meshWidth, parti)
{};


double CentralDifferences::computeDu2Dx (int i, int j) const
{
	return ((u_(i+1,j)+u_(i,j))*(u_(i+1,j)+u_(i,j))-(u_(i,j)+u_(i-1,j))*(u_(i,j)+u_(i-1,j)))/(4*meshWidth_[0]);
}

double CentralDifferences::computeDv2Dy (int i, int j) const
{
	return ((v_(i,j+1)+v_(i,j))*(v_(i,j+1)+v_(i,j))-(v_(i,j)+v_(i,j-1))*(v_(i,j)+v_(i,j-1)))/(4*meshWidth_[1]);
}

double CentralDifferences::computeDuvDx (int i, int j) const
{
	return ((u_(i,j+1)+u_(i,j))*(v_(i,j)+v_(i+1,j))-(u_(i-1,j+1)+u_(i-1,j))*(v_(i-1,j)+v_(i,j)))/(4*meshWidth_[0]);
}

double CentralDifferences::computeDuvDy (int i, int j) const
{
	return ((v_(i+1,j)+v_(i,j))*(u_(i,j+1)+u_(i,j))-(v_(i+1,j-1)+v_(i,j-1))*(u_(i,j)+u_(i,j-1)))/(4*meshWidth_[1]);
}
