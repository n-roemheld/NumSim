#include "CentralDifferences.h"


// CentralDifferences::CentralDifferences(std::array< int, 2 > nCells, std::array< double, 2 > meshWidth) :
// 	Discretization(nCells, meshWidth)
// {};

// CentralDifferences::CentralDifferences(std::array< int, 2 > nCells, std::array< double, 2 > meshWidth, std::shared_ptr<Array2D> geometryPVString, std::shared_ptr<Array2D> geometryPV1, std::shared_ptr<Array2D> geometryPV2, std::shared_ptr<Array2D> geometryTString, std::shared_ptr<Array2D> geometryT1) :
// 	Discretization(nCells, meshWidth, geometryPVString, geometryPV1, geometryPV2, geometryTString, geometryT1)
// {};

CentralDifferences::CentralDifferences(std::array< int, 2 > nCells, std::array< double, 2 > meshWidth, std::shared_ptr<Array2D> geometryPVString, std::shared_ptr<Array2D> geometryPV1, std::shared_ptr<Array2D> geometryPV2, std::shared_ptr<Array2D> geometryTString, std::shared_ptr<Array2D> geometryT1, Adapter& adapter) :
	Discretization(nCells, meshWidth, geometryPVString, geometryPV1, geometryPV2, geometryTString, geometryT1, adapter)
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

double CentralDifferences::computeDuTDx(int i, int j) const
{
	return (u(i,j)*(T(i+1,j)+T(i,j))/2-u(i-1,j)*(T(i,j)+T(i-1,j))/2)/meshWidth_[0];
}

double CentralDifferences::computeDvTDy(int i, int j) const
{
	return (v(i,j)*(T(i,j+1)+T(i,j))/2-v(i,j-1)*(T(i,j)+T(i,j-1))/2)/meshWidth_[1];
}
