#include "DonorCell.h"
#include <math.h>

// DonorCell::DonorCell(std::array<int, 2> nCells, std::array<double, 2> meshWidth,
// 		double alpha) :
// 		Discretization(nCells, meshWidth), alpha_(alpha)
// {};

// DonorCell::DonorCell(std::array< int, 2 > nCells, std::array< double, 2 > meshWidth, std::shared_ptr<Array2D> geometryPVString, std::shared_ptr<Array2D> geometryPV1, std::shared_ptr<Array2D> geometryPV2, std::shared_ptr<Array2D> geometryTString, std::shared_ptr<Array2D> geometryT1, double alpha, double gamma) :
// 	Discretization(nCells, meshWidth, geometryPVString, geometryPV1, geometryPV2, geometryTString, geometryT1), alpha_(alpha), gamma_(gamma)
// {};

DonorCell::DonorCell(std::array< int, 2 > nCells, std::array< double, 2 > meshWidth, std::shared_ptr<Array2D> geometryPVString, std::shared_ptr<Array2D> geometryPV1, std::shared_ptr<Array2D> geometryPV2, std::shared_ptr<Array2D> geometryTString, std::shared_ptr<Array2D> geometryT1, double alpha, double gamma, Adapter& adapter) :
	Discretization(nCells, meshWidth, geometryPVString, geometryPV1, geometryPV2, geometryTString, geometryT1, adapter), alpha_(alpha), gamma_(gamma)
{};

// alpha parameterized donor cell /centrall difference scheme
double DonorCell::computeDu2Dx(int i, int j) const {

    return 1 / dx() * (pow((u(i, j) + u(i + 1, j)) / 2, 2) - pow((u(i - 1, j) + u(i, j)) / 2, 2))
           + alpha_ * 1 / dx() * (fabs(u(i, j) + u(i + 1, j)) / 2 * (u(i, j) - u(i + 1, j)) / 2 -
                                  fabs(u(i - 1, j) + u(i, j)) / 2 * (u(i - 1, j) - u(i, j)) / 2);
}

// alpha parameterized donor cell /centrall difference scheme
double DonorCell::computeDv2Dy(int i, int j) const {
    return 1 / dy() * (pow((v(i, j) + v(i, j + 1)) / 2, 2) - pow((v(i, j - 1) + v(i, j)) / 2, 2))
           + alpha_ * 1 / dy() * (fabs(v(i, j) + v(i, j + 1)) / 2 * (v(i, j) - v(i, j + 1)) / 2 -
                                  fabs(v(i, j - 1) + v(i, j)) / 2 * (v(i, j - 1) - v(i, j)) / 2);
}

// alpha parameterized donor cell /centrall difference scheme
double DonorCell::computeDuvDx(int i, int j) const {
    return 1 / dx() * ((u(i, j) + u(i, j + 1)) / 2 * (v(i, j) + v(i + 1, j)) / 2 -
                       (u(i - 1, j) + u(i - 1, j + 1)) / 2 * (v(i - 1, j) + v(i, j)) / 2)
           + alpha_ * 1 / dx() * (fabs(u(i, j) + u(i, j + 1)) / 2 * (v(i, j) - v(i + 1, j)) / 2 -
                                  fabs(u(i - 1, j) + u(i - 1, j + 1)) / 2 * (v(i - 1, j) - v(i, j)) / 2);
}

// alpha parameterized donor cell /centrall difference scheme
double DonorCell::computeDuvDy(int i, int j) const {
    return 1 / dy() * ((v(i, j) + v(i + 1, j)) / 2 * (u(i, j) + u(i, j + 1)) / 2 -
                       (v(i, j - 1) + v(i + 1, j - 1)) / 2 * (u(i, j - 1) + u(i, j)) / 2)
           + alpha_ * 1 / dy() * (fabs(v(i, j) + v(i + 1, j)) / 2 * (u(i, j) - u(i, j + 1)) / 2 -
                                  fabs(v(i, j - 1) + v(i + 1, j - 1)) / 2 * (u(i, j - 1) - u(i, j)) / 2);
}

double DonorCell::computeDuTDx(int i, int j) const
{
	return 1/meshWidth_[0]*(u(i,j)*(T(i+1,j)+T(i,j))/2-u(i-1,j)*(T(i,j)+T(i-1,j))/2
                         + gamma_*( fabs(u(i,j))*  (T(i,j)-T(i+1,j))/2
			                         -fabs(u(i-1,j))*(T(i-1,j)-T(i,j))/2 ));
}

double DonorCell::computeDvTDy(int i, int j) const
{
	return 1/meshWidth_[1]*(v(i,j)*(T(i,j+1)+T(i,j))/2-v(i,j-1)*(T(i,j)+T(i,j-1))/2
                         + gamma_*( fabs(v(i,j))  *(T(i,j)-T(i,j+1))/2
			                        - fabs(v(i,j-1))*(T(i,j-1)-T(i,j))/2 ));
}
