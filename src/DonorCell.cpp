#include "DonorCell.h"
#include <math.h>

DonorCell::DonorCell(std::array<int, 2> nCells, std::array<double, 2> meshWidth,
		double alpha) :
		Discretization(nCells, meshWidth), alpha_(alpha)
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
//    return 1 / dx() * ((u(i, j) + u(i, j + 1)) / 2 * (v(i, j) + v(i + 1, j)) / 2 -
//                       (u(i - 1, j) + u(i - 1, j + 1)) / 2 * (v(i - 1, j) + v(i, j)) / 2)
//           + alpha_ * 1 / dx() * (fabs(u(i, j) + u(i, j + 1)) / 2 * (v(i, j) - v(i + 1, j)) / 2 -
//                                  fabs(u(i - 1, j) + u(i - 1, j + 1)) / 2 * (v(i - 1, j) - v(i, j)) / 2);
	   return 1 / dx() * ((u(i+1, j) + u(i+1, j + 1)) / 2 * (v(i, j) + v(i + 1, j)) / 2 -
	                       (u(i, j) + u(i, j + 1)) / 2 * (v(i - 1, j) + v(i, j)) / 2)
	           + alpha_ * 1 / dx() * (fabs(u(i+1, j) + u(i+1, j + 1)) / 2 * (v(i, j) - v(i + 1, j)) / 2 -
	                                  fabs(u(i, j) + u(i , j + 1)) / 2 * (v(i - 1, j) - v(i, j)) / 2);
}

// alpha parameterized donor cell /centrall difference scheme
double DonorCell::computeDuvDy(int i, int j) const {
//    return 1 / dy() * ((v(i, j) + v(i + 1, j)) / 2 * (u(i, j) + u(i, j + 1)) / 2 -
//                       (v(i, j - 1) + v(i + 1, j - 1)) / 2 * (u(i, j - 1) + u(i, j)) / 2)
//           + alpha_ * 1 / dy() * (fabs(v(i, j) + v(i + 1, j)) / 2 * (u(i, j) - u(i, j + 1)) / 2 -
//                                  fabs(v(i, j - 1) + v(i + 1, j - 1)) / 2 * (u(i, j - 1) - u(i, j)) / 2);
	return 1 / dy() * ((v(i, j+1) + v(i + 1, j+1)) / 2 * (u(i, j) + u(i, j + 1)) / 2 -
	                       (v(i, j) + v(i + 1, j)) / 2 * (u(i, j - 1) + u(i, j)) / 2)
	           + alpha_ * 1 / dy() * (fabs(v(i, j+1) + v(i + 1, j+1)) / 2 * (u(i, j) - u(i, j + 1)) / 2 -
	                                  fabs(v(i, j) + v(i + 1, j)) / 2 * (u(i, j - 1) - u(i, j)) / 2);
}
