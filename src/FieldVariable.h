#pragma once

#include <array>
#include "Array2D.h"

class FieldVariable : public Array2D
{
public:

	//!constructor
	FieldVariable (std::array< int, 2 > size, std::array< double, 2 > origin, std::array< double, 2 > meshWidth);

	//!get the value at the Cartesian coordinate (x,y). The value is linearly interpolated between stored points.
	double interpolateAt(double x, double y) const;

	FieldVariable& operator=(FieldVariable rhs);

private:
	// to do: make const again!
	const std::array<double,2> origin_;
	// to do: make const again!
	const std::array<double,2> meshWidth_;

};
