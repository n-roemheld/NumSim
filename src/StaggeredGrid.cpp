#include "StaggeredGrid.h"

StaggeredGrid::StaggeredGrid(std::array<int,2> nCells, std::array<double,2> meshWidth) :
	nCells_(nCells), meshWidth_(meshWidth),
	u_( {nCells[0]+1, nCells[1]+2},  {0*meshWidth[0],    -0.5*meshWidth[1]}, meshWidth ),
	v_( {nCells[0]+2, nCells[1]+1},  {-0.5*meshWidth[0], 0*meshWidth[1]}, meshWidth ),
	p_( {nCells[0]+2, nCells[1]+2},  {-0.5*meshWidth[0], -0.5*meshWidth[1]}, meshWidth),
	f_( {nCells[0]+1, nCells[1]+2},  {0*meshWidth[0],    -0.5*meshWidth[1]}, meshWidth ),
	g_( {nCells[0]+2, nCells[1]+1},  {-0.5*meshWidth[0], 0*meshWidth[1]}, meshWidth ),
	rhs_( {nCells[0]+2, nCells[1]+2},{-0.5*meshWidth[0], -0.5*meshWidth[1]}, meshWidth),
	T_( {nCells[0]+2, nCells[1]+2},  {-0.5*meshWidth[0], -0.5*meshWidth[1]}, meshWidth)
{};

StaggeredGrid::StaggeredGrid(std::array< int, 2 > nCells, std::array< double, 2 > meshWidth, std::shared_ptr<Array2D> geometryPVString, std::shared_ptr<Array2D> geometryPV1, std::shared_ptr<Array2D> geometryPV2, std::shared_ptr<Array2D> geometryTString, std::shared_ptr<Array2D> geometryT1) :
	nCells_(nCells), meshWidth_(meshWidth),
	u_( {nCells[0]+1, nCells[1]+2},  {0*meshWidth[0],    -0.5*meshWidth[1]}, meshWidth ),
	v_( {nCells[0]+2, nCells[1]+1},  {-0.5*meshWidth[0], 0*meshWidth[1]}, meshWidth ),
	p_( {nCells[0]+2, nCells[1]+2},  {-0.5*meshWidth[0], -0.5*meshWidth[1]}, meshWidth),
	f_( {nCells[0]+1, nCells[1]+2},  {0*meshWidth[0],    -0.5*meshWidth[1]}, meshWidth ),
	g_( {nCells[0]+2, nCells[1]+1},  {-0.5*meshWidth[0], 0*meshWidth[1]}, meshWidth ),
	rhs_( {nCells[0]+2, nCells[1]+2},{-0.5*meshWidth[0], -0.5*meshWidth[1]}, meshWidth),
	T_( {nCells[0]+2, nCells[1]+2},  {-0.5*meshWidth[0], -0.5*meshWidth[1]}, meshWidth),
	geometryPVString_(geometryPVString), geometryPV1_(geometryPV1), geometryPV2_(geometryPV2), geometryTString_(geometryTString), geometryT1_(geometryT1)
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

// removed const after function
const FieldVariable& StaggeredGrid::T() const
{
  return T_;
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

double& StaggeredGrid::T(int i, int j)
{
  return StaggeredGrid::T_(i,j);
};

double StaggeredGrid::T(int i, int j) const
{
  return StaggeredGrid::T_(i,j);
};

double StaggeredGrid::geometryPVString(int i, int j) const
{
	return StaggeredGrid::geometryPVString_->operator()(i, j);
}

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
	return 1;
};

int StaggeredGrid::uIEnd() const
{
	return nCells_[0]+1; // reduced by one to account for boundary values // or not? // or not not?
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
	return 1;
};

int StaggeredGrid::vJEnd() const
{
	return nCells_[1]+1; // reduced by one to account for boundary values // or not? // or not not?
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


void StaggeredGrid::setBoundaryValues_u_f(int location_boundary, int i, int j)
{
	// indices in geometry file (shifted by uIBegin and increased by one at the right (u) and upper(v) boundaries)
	int igeom = i-uIBegin()+1; // todo: double check!!
	int jgeom = j-uJBegin()+1;
	// indices of the neigboring cell (can be fluid or solid)
	int in = i;
	int jn = j;
	switch (location_boundary) // 0: Left, 1: Right, 2: Lower, 3: Upper
	{
		case 0:	in = i+1; igeom = 0; break;
		case 1: in = i-1; igeom = geometryPVString_->size()[0]-1; break;
		case 2: jn = j+1; jgeom = 0; break;
		case 3: jn = j-1; jgeom = geometryPVString_->size()[1]-1; break;
	}
	// set boundary values to nan if not needed (neighbor not fluid cell)
	if (false){//((geometryPVString_->operator()(igeom + in - i, jgeom + jn - j) != -1)	{
		// std::cout << "nan1: " << geometryPVString_->operator()(igeom + in - i, jgeom + jn - j) << " " << igeom + in - i << " " << jgeom + jn - j << " " << igeom << i << in << " " << jgeom << j << jn << std::endl;
		// u(i,j) = std::nan("1"); commented out by Henrik
	}
	else
	{
		// set boundary values of one cell
		int boundary_type = geometryPVString_->operator()(igeom, jgeom);
		switch (boundary_type)
		{
			case 0: // NOSLIP
			{
				if (location_boundary == 0 || location_boundary == 1) { // left or right
					u(i,j) = 0;
				} else { // upper or lower
					u(i,j) = - u(in,jn);
				}
				break;
			}
			case 1: // SLIP
			{
				if (location_boundary == 0 || location_boundary == 1) { // left or right
					u(i,j) = 0;
				} else { // upper or lower
					u(i,j) = u(in,jn);
				}
				break;
			}
			case 2: // INFLOW
			{
				if (location_boundary == 0 || location_boundary == 1) { // left or right
					u(i,j) = geometryPV1_->operator()(igeom,jgeom);
				} else { // upper or lower
					u(i,j) = 2*geometryPV1_->operator()(igeom,jgeom) - u(in,jn);
				}
				break;
			}
			case 3: // OUTFLOW (same as case 4)
			case 4: // PRESSURE
			{
				if (location_boundary == 0 || location_boundary == 1) { // left or right
					u(i,j) = u(in,jn);
				} else { // upper or lower
					u(i,j) = u(in,jn);
				}
				break;
			}

		}
	}
	f(i,j) = u(i,j);
}

void StaggeredGrid::setBoundaryValues_v_g(int location_boundary, int i, int j)
{
	// indices in geometry file (shifted by uIBegin and increased by one at the right (u) and upper(v) boundaries)
	int igeom = i-vIBegin()+1; // todo: double check!!
	int jgeom = j-vJBegin()+1;
	// indices of the neigboring cell (can be fluid or solid)
	int in = i;
	int jn = j;
	switch (location_boundary) // 0: Left, 1: Right, 2: Lower, 3: Upper
	{
		case 0:	in = i+1; igeom = 0; break;
		case 1: in = i-1; igeom = geometryPVString_->size()[0]-1; break;
		case 2: jn = j+1; jgeom = 0; break;
		case 3: jn = j-1; jgeom = geometryPVString_->size()[1]-1; break;
	}
	// set boundary values to nan if not needed (neighbor not fluid cell)
	if (false){//(geometryPVString_->operator()(igeom + in - i, jgeom + jn - j) != -1)	{
		// v(i,j) = std::nan("1"); commented out by Henrik
	}
	else
	{
		// set boundary values of one cell
		int boundary_type = geometryPVString_->operator()(igeom, jgeom);
		switch (boundary_type)
		{
			case 0: // NOSLIP
			{
				if (location_boundary == 0 || location_boundary == 1) { // left or right
					v(i,j) = - v(in,jn);
				} else { // upper or lower
					v(i,j) = 0;
				}
				break;
			}
			case 1: // SLIP
			{
				if (location_boundary == 0 || location_boundary == 1) { // left or right
					v(i,j) = v(in,jn);
				} else { // upper or lower
					v(i,j) = 0;
				}
				break;
			}
			case 2: // INFLOW
			{
				if (location_boundary == 0 || location_boundary == 1) { // left or right
					v(i,j) = 2*geometryPV2_->operator()(igeom,jgeom) - v(in,jn);
				} else { // upper or lower
					v(i,j) = geometryPV2_->operator()(igeom,jgeom);
				}
				break;
			}
			case 3: // OUTFLOW (same as case 4)
			case 4: // PRESSURE
			{
				if (location_boundary == 0 || location_boundary == 1) { // left or right
					v(i,j) = v(in,jn);
				} else { // upper or lower
					v(i,j) = v(in,jn);
				}
				break;
			}

		}
	}
	g(i,j) = v(i,j);
}

void StaggeredGrid::setBoundaryValues_p(int location_boundary, int i, int j)
{
	// indices in geometry file (shifted by uIBegin and increased by one at the right (u) and upper(v) boundaries)
	int igeom = i-pIBegin()+1; // todo: double check!!
	int jgeom = j-pJBegin()+1;
	// indices of the neigboring cell (can be fluid or solid)
	int in = i;
	int jn = j;
	switch (location_boundary) // 0: Left, 1: Right, 2: Lower, 3: Upper
	{
		case 0:	in = i+1; igeom = 0; break;
		case 1: in = i-1; igeom = geometryPVString_->size()[0]-1; break;
		case 2: jn = j+1; jgeom = 0; break;
		case 3: jn = j-1; jgeom = geometryPVString_->size()[1]-1; break;
	}
	// set boundary values to nan if not needed (neighbor not fluid cell)
	if (false){//((geometryPVString_->operator()(igeom + in - i, jgeom + jn - j) != -1)	{
		// std::cout << "nan3: " << location_boundary << ", " << in << ", " << i << ", " << igeom << ", j: " << jn << j << jgeom << std::endl;
		// p(i,j) = std::nan("1"); commented out by Henrik
	}
	else
	{
		// set boundary values of one cell
		int boundary_type = geometryPVString_->operator()(igeom, jgeom);
		if (boundary_type == 4) // PRESSURE
		{
			p(i,j) = 2*geometryPV1_->operator()(igeom, jgeom) - p(in,jn);
		}
		else
		{
			p(i,j) = p(in,jn);
		}

	}
}

void StaggeredGrid::setBoundaryValues_T(int location_boundary, int i, int j)
{
	// indices in geometry file (shifted by uIBegin and increased by one at the right (u) and upper(v) boundaries)
	int igeom = i-pIBegin()+1; // todo: double check!!
	int jgeom = j-pJBegin()+1;
	// mesh width h
	double h = meshWidth_[0];
	// indices of the neigboring cell (can be fluid or solid)
	int in = i;
	int jn = j;
	switch (location_boundary) // 0: Left, 1: Right, 2: Lower, 3: Upper
	{
		case 0:	in = i+1; igeom = 0; break;
		case 1: in = i-1; igeom = geometryPVString_->size()[0]-1; break;
		case 2: jn = j+1; jgeom = 0; h = meshWidth_[1]; break;
		case 3: jn = j-1; jgeom = geometryPVString_->size()[1]-1; h = meshWidth_[1]; break;
	}
	// set boundary values to nan if not needed (neighbor not fluid cell)
	if (false){//((geometryPVString_->operator()(igeom + in - i, jgeom + jn - j) != -1)	{
		// T(i,j) = std::nan("1"); commented out by Henrik
	}
	else
	{
		// set boundary values of one cell
		int boundary_type = geometryTString_->operator()(igeom, jgeom);
		if (boundary_type == 0) // TD
		{
			T(i,j) = 2*geometryT1_->operator()(igeom, jgeom) - T(in,jn);
		}
		else // TN
		{
			T(i,j) = T(in,jn) + h*geometryT1_->operator()(igeom, jgeom);
		}

	}
}

void StaggeredGrid::fillIn(int uInit, int vInit, int pInit, int TInit)
{
	  for(int j = 0; j < u_.size()[1]; j++)
	  {
	    for(int i = 0; i < u_.size()[0]; i++)
	    {
	      u_(i,j) = uInit;
	    }
	  }

		for(int j = 0; j < v_.size()[1]; j++)
		{
			for(int i = 0; i < v_.size()[0]; i++)
			{
				v_(i,j) = vInit;
			}
		}

		for(int j = 0; j < p_.size()[1]; j++)
		{
			for(int i = 0; i < p_.size()[0]; i++)
			{
				p_(i,j) = pInit;
				T_(i,j) = TInit;
			}
		}
}
