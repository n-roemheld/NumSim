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

// set obstaclesValues iff cell is solid for u and f
void StaggeredGrid::setObstacleValues_u_f(int i, int j)
{
	std::array<int , 2> locations_boundary = {-2,-2};
	// indices in geometry file (shifted by uIBegin and increased by one at the right (u) and upper(v) boundaries)
	int igeom = i-uIBegin()+1; // todo: double check!!
	int jgeom = j-uJBegin()+1;
	if(geometryPVString_->operator()(igeom, jgeom) == 5) // proof, if cell is solid
	{
	int index = 0;
	for(int i = 0; i < 4; i++)
	{
		if(index > 1)
		{
			std::cout << "u: boundary locations are not calculated correctly or two-cell-criterion is not fullfilled" << std::endl;
		}
		
		switch(i) // 0: Left, 1: Right, 2: Lower, 3: Upper
		{
			case 0:
				if(geometryPVString_->operator()(igeom-1, jgeom) == -1)
				{
					locations_boundary[index] = i;
					index++;
				}
				break;
			case 1:
				if(geometryPVString_->operator()(igeom+1, jgeom) == -1)
				{
					locations_boundary[index] = i;
					index++;
				}
				break;
			case 2:
				if(geometryPVString_->operator()(igeom, jgeom-1) == -1)
				{
					locations_boundary[index] = i;
					index++;
				}
				break;
			case 3:
				if(geometryPVString_->operator()(igeom, jgeom+1) == -1)
				{
					locations_boundary[index] = i;
					index++;
				}
				break;
		}
	}

	while(index < 2)
	{
		locations_boundary[index] = -1;
		index++;
	}

	setObstacleValues_u_f(locations_boundary, i, j);
	}
	
}

void StaggeredGrid::setObstacleValues_v_g(int i, int j)
{
	std::array<int , 2> locations_boundary = {-2,-2};
	// indices in geometry file (shifted by uIBegin and increased by one at the right (u) and upper(v) boundaries)
	int igeom = i-vIBegin()+1; // todo: double check!!
	int jgeom = j-vJBegin()+1;
	if(geometryPVString_->operator()(igeom, jgeom) == 5) // proof, if cell is solid
	{
	int index = 0;
	for(int i = 0; i < 4; i++)
	{
		if(index > 1)
		{
			std::cout << "v: boundary locations are not calculated correctly or two-cell-criterion is not fullfilled" << std::endl;
		}
		switch(i) // 0: Left, 1: Right, 2: Lower, 3: Upper
		{
			case 0:
				if(geometryPVString_->operator()(igeom-1, jgeom) == -1)
				{
					locations_boundary[index] = i;
					index++;
				}
				break;
			case 1:
				if(geometryPVString_->operator()(igeom+1, jgeom) == -1)
				{
					locations_boundary[index] = i;
					index++;
				}
				break;
			case 2:
				if(geometryPVString_->operator()(igeom, jgeom-1) == -1)
				{
					locations_boundary[index] = i;
					index++;
				}
				break;
			case 3:
				if(geometryPVString_->operator()(igeom, jgeom+1) == -1)
				{
					locations_boundary[index] = i;
					index++;
				}
				break;
		}
	}

	while(index < 2)
	{
		locations_boundary[index] = -1;
		index++;
	}

	setObstacleValues_v_g(locations_boundary, i, j);
	}
}

void StaggeredGrid::setObstacleValues_p(int i, int j)
{
	std::array<int , 2> locations_boundary = {-2,-2};
	// indices in geometry file (shifted by uIBegin and increased by one at the right (u) and upper(v) boundaries)
	int igeom = i-pIBegin()+1; // todo: double check!!
	int jgeom = j-pJBegin()+1;
	if(geometryPVString_->operator()(igeom, jgeom) == 5) // proof, if cell is solid
	{
	int index = 0;
	for(int i = 0; i < 4; i++)
	{
		if(index > 1)
		{
			std::cout << "p: boundary locations are not calculated correctly or two-cell-criterion is not fullfilled" << std::endl;
		}
		switch(i) // 0: Left, 1: Right, 2: Lower, 3: Upper
		{
			case 0:
				if(geometryPVString_->operator()(igeom-1, jgeom) == -1)
				{
					locations_boundary[index] = i;
					index++;
				}
				break;
			case 1:
				if(geometryPVString_->operator()(igeom+1, jgeom) == -1)
				{
					locations_boundary[index] = i;
					index++;
				}
				break;
			case 2:
				if(geometryPVString_->operator()(igeom, jgeom-1) == -1)
				{
					locations_boundary[index] = i;
					index++;
				}
				break;
			case 3:
				if(geometryPVString_->operator()(igeom, jgeom+1) == -1)
				{
					locations_boundary[index] = i;
					index++;
				}
				break;
		}
	}

	while(index < 2)
	{
		locations_boundary[index] = -1;
		index++;
	}

	setObstacleValues_p(locations_boundary, i, j);
	}
}

void StaggeredGrid::setObstacleValues_T(int i, int j)
{
	std::array<int , 2> locations_boundary = {-2,-2};
	// indices in geometry file (shifted by uIBegin and increased by one at the right (u) and upper(v) boundaries)
	int igeom = i-pIBegin()+1; // todo: double check!!
	int jgeom = j-pJBegin()+1;
	if(geometryPVString_->operator()(igeom, jgeom) == 5) // proof, if cell is solid
	{
	int index = 0;
	for(int i = 0; i < 4; i++)
	{
		if(index > 1)
		{
			std::cout << "T: boundary locations are not calculated correctly or two-cell-criterion is not fullfilled" << std::endl;
		}
		switch(i) // 0: Left, 1: Right, 2: Lower, 3: Upper
		{
			case 0:
				if(geometryPVString_->operator()(igeom-1, jgeom) == -1)
				{
					locations_boundary[index] = i;
					index++;
				}
				break;
			case 1:
				if(geometryPVString_->operator()(igeom+1, jgeom) == -1)
				{
					locations_boundary[index] = i;
					index++;
				}
				break;
			case 2:
				if(geometryPVString_->operator()(igeom, jgeom-1) == -1)
				{
					locations_boundary[index] = i;
					index++;
				}
				break;
			case 3:
				if(geometryPVString_->operator()(igeom, jgeom+1) == -1)
				{
					locations_boundary[index] = i;
					index++;
				}
				break;
		}
	}

	while(index < 2)
	{
		locations_boundary[index] = -1;
		index++;
	}

	setObstacleValues_T(locations_boundary, i, j);
	}
}

void StaggeredGrid::setObstacleValues_u_f2(int i, int j)
{
	int igeom = i-uIBegin()+1; // todo: double check!!
	int jgeom = j-uJBegin()+1;

	if (geometryPVString_->operator()(igeom - 1, jgeom) == -1) // left is fluid
	{
		u_(i-1,j) = 0; f_(i-1,j) = u_(i-1,j);
		if (geometryPVString_->operator()(igeom, jgeom - 1) == -1) // left and lower is fluid //todo: correct grammar
		{
			u_(i,j) = -u_(i,j-1);
		}
		else if (geometryPVString_->operator()(igeom, jgeom + 1) == -1) // left and upper is fluid
		{
			u_(i,j) = -u_(i,j+1);
		}
		else // only left is fluid
		{
			u_(i,j) = std::nan("1");
		}
	}
	else if (geometryPVString_->operator()(igeom + 1, jgeom) == -1) // right is fluid
	{
		u_(i,j) = 0;
	}
	else // left and right not fluid
	{
		if (geometryPVString_->operator()(igeom, jgeom - 1) == -1) // only lower is fluid
		{
			u_(i,j) = -u_(i,j-1);
		}
		else if (geometryPVString_->operator()(igeom, jgeom + 1) == -1) // only upper is fluid
		{
			u_(i,j) = -u_(i,j+1);
		}
		else // no direct neighbors are fluid
		{
			u_(i,j) = std::nan("1");
		}
	}

	f_(i,j) = u(i,j);

	if (geometryPVString_->operator()(igeom, jgeom + 1) == -1   &&   geometryPVString_->operator()(igeom + 1, jgeom + 1) != -1) // upper is fluid but not upper-right
	{
		u_(i,j) = std::nan("1");
	}
	else if (geometryPVString_->operator()(igeom, jgeom - 1) == -1   &&   geometryPVString_->operator()(igeom + 1, jgeom - 1) != -1) // lower is fluid, but not lower-right
	{
		u_(i,j) = std::nan("1");
	}
	else if (geometryPVString_->operator()(igeom - 1, jgeom) == -1   &&   geometryPVString_->operator()(igeom - 1, jgeom + 1) != -1) // left is fluid, but not upper-left
	{
		u_(i,j) = std::nan("1");
	}
}

void StaggeredGrid::setObstacleValues_v_g2(int i, int j)
{
	int igeom = i-vIBegin()+1; // todo: double check!!
	int jgeom = j-vJBegin()+1;

	if (geometryPVString_->operator()(igeom - 1, jgeom) == -1) // left is fluid
	{
		if (geometryPVString_->operator()(igeom, jgeom - 1) == -1) // left and lower is fluid //todo: correct grammar
		{
			v_(i,j-1) = 0; g_(i,j-1) = v_(i,j-1);
			v_(i,j) = -v_(i-1,j);
		}
		else if (geometryPVString_->operator()(igeom, jgeom + 1) == -1) // left and upper is fluid
		{
			v_(i,j) = 0;
		}
		else // only left is fluid
		{
			v_(i,j) = -v_(i-1,j);
		}
	}
	else if (geometryPVString_->operator()(igeom + 1, jgeom) == -1) // right is fluid
	{
		v_(i,j) = -v_(i+1,j);
		if (geometryPVString_->operator()(igeom, jgeom - 1) == -1) // right and lower is fluid
		{
			v_(i,j-1) = 0; g_(i,j-1) = v_(i,j-1);
		}	
	}
	else // left and right not fluid
	{
		if (geometryPVString_->operator()(igeom, jgeom - 1) == -1) // only lower is fluid
		{
			v_(i,j-1) = 0; g_(i,j-1) = v_(i,j-1);
			v_(i,j) = std::nan("1");
		}
		else if (geometryPVString_->operator()(igeom, jgeom + 1) == -1) // only upper is fluid
		{
			v_(i,j) = 0;
		}
		else // no direct neighbors are fluid
		{
			v_(i,j) = std::nan("1");
		}
	}

	g_(i,j) = v(i,j);
	if (geometryPVString_->operator()(igeom, jgeom - 1) == -1   &&   geometryPVString_->operator()(igeom + 1, jgeom - 1) != -1) // lower is fluid, but not lower-right
	{
		v_(i,j) = std::nan("1");
	}
	else if (geometryPVString_->operator()(igeom + 1, jgeom) == -1   &&   geometryPVString_->operator()(igeom + 1, jgeom + 1) != -1) // right is fluid but not upper-right
	{
		v_(i,j) = std::nan("1");
	}
	else if (geometryPVString_->operator()(igeom - 1, jgeom) == -1   &&   geometryPVString_->operator()(igeom - 1, jgeom + 1) != -1) // left is fluid, but not upper-left
	{
		v_(i,j) = std::nan("1");
	}
}

void StaggeredGrid::setObstacleValues_p2(int i, int j)
{
	int igeom = i-pIBegin()+1; // todo: double check!!
	int jgeom = j-pJBegin()+1;

	if (geometryPVString_->operator()(igeom, jgeom) == 5) // if current cell is not solid
	{
		return;
	}

	if (geometryPVString_->operator()(igeom - 1, jgeom) == -1) // left is fluid
	{
		if (geometryPVString_->operator()(igeom, jgeom - 1) == -1) // left and lower is fluid //todo: correct grammar
		{
			p_(i,j) = .5*(p_(i-1,j) + p_(i,j-1));
		}
		else if (geometryPVString_->operator()(igeom, jgeom + 1) == -1) // left and upper is fluid
		{
		}
		else // only left is fluid
		{
			p_(i,j) = p_(i-1,j);
		}
	}
	else if (geometryPVString_->operator()(igeom + 1, jgeom) == -1) // right is fluid
	{
		if (geometryPVString_->operator()(igeom, jgeom - 1) == -1) // right and lower is fluid
		{
			p_(i,j) = .5*(p_(i+1,j) + p_(i,j-1));
		}
		else if (geometryPVString_->operator()(igeom, jgeom + 1) == -1) // right and upper is fluid
		{
			p_(i,j) = .5*(p_(i+1,j) + p_(i,j+1));
		}
		else // only right is fluid
		{
			p_(i,j) = p_(i+1,j);
		}		
	}
	else // left and right not fluid
	{
		if (geometryPVString_->operator()(igeom, jgeom - 1) == -1) // only lower is fluid
		{
			p_(i,j) = p_(i,j-1);
		}
		else if (geometryPVString_->operator()(igeom, jgeom + 1) == -1) // only upper is fluid
		{
			p_(i,j) = p_(i,j+1);
		}
		else // no direct neighbors are fluid
		{
			p_(i,j) = std::nan("1");
		}
		
	}
}

void StaggeredGrid::setObstacleValues_T2(int i, int j)
{
	int igeom = i-pIBegin()+1; // todo: double check!!
	int jgeom = j-pJBegin()+1;

	if (geometryPVString_->operator()(igeom - 1, jgeom) == -1) // left is fluid
	{
		if (geometryPVString_->operator()(igeom, jgeom - 1) == -1) // left and lower is fluid //todo: correct grammar
		{
			T_(i,j) = .5*(T_(i-1,j) + T_(i,j-1));
		}
		else if (geometryPVString_->operator()(igeom, jgeom + 1) == -1) // left and upper is fluid
		{
		}
		else // only left is fluid
		{
			T_(i,j) = T_(i-1,j);
		}
	}
	else if (geometryPVString_->operator()(igeom + 1, jgeom) == -1) // right is fluid
	{
		if (geometryPVString_->operator()(igeom, jgeom - 1) == -1) // right and lower is fluid
		{
			T_(i,j) = .5*(T_(i+1,j) + T_(i,j-1));
		}
		else if (geometryPVString_->operator()(igeom, jgeom + 1) == -1) // right and upper is fluid
		{
			T_(i,j) = .5*(T_(i+1,j) + T_(i,j+1));
		}
		else // only right is fluid
		{
			T_(i,j) = T_(i+1,j);
		}		
	}
	else // left and right not fluid
	{
		if (geometryPVString_->operator()(igeom, jgeom - 1) == -1) // only lower is fluid
		{
			T_(i,j) = T_(i,j-1);
		}
		else if (geometryPVString_->operator()(igeom, jgeom + 1) == -1) // only upper is fluid
		{
			T_(i,j) = T_(i,j+1);
		}
		else // no direct neighbors are fluid
		{
			T_(i,j) = std::nan("1");
		}
		
	}
}


void StaggeredGrid::setObstacleValues_u_f(std::array<int, 2> locations_boundary, int i, int j)
{
	// sort locations_boundary from big to small
	if(locations_boundary[0] < locations_boundary[1])
	{
		int temp = locations_boundary[0];
		locations_boundary[0] = locations_boundary[1];
		locations_boundary[1] = temp;
	}

	// set value to nan if cell is inner obstacal cell
	if(locations_boundary[0] == -1 && locations_boundary[1] == -1)
	{
		u_(i,j) = std::nan("1");
		f_(i,j) = u_(i,j);
	}

	// proof if corner or not 
	if(locations_boundary[1] == -1)
	{
		switch(locations_boundary[0]) // 0: Left, 1: Right, 2: Lower, 3: Upper
		{
			case 0:
				u_(i-1,j) = 0;
				f_(i-1,j) = u_(i-1,j);
				u_(i,j) = std::nan("1");
				f_(i,j) = u_(i,j);
				break;
			case 1:
				u_(i,j) = 0;
				f_(i,j) = u_(i,j);
				break;
			case 2:
				u_(i,j) = -u_(i,j-1);
				f_(i,j) = u_(i,j);
				break;
			case 3:
				u_(i,j) = -u_(i,j+1);
				f_(i,j) = u_(i,j);
				break;
		}
	}
	else 
	{
		if(locations_boundary[0] == 2)
		{
			switch(locations_boundary[1])
			{
				case 0: //edge lower left
					u_(i-1,j) = 0;
					f_(i-1,j) = u_(i-1,j);
					u_(i,j) = -u_(i,j-1);
					f_(i,j) = u_(i,j);
					break;
				case 1: //edge lower right
					u_(i,j) = 0;
					f_(i,j) = u_(i,j);
					break;
			}
		} 
		else if(locations_boundary[0] == 3)
		{
			switch(locations_boundary[1])
			{
				case 0: //edge upper left
					u_(i-1,j) = 0;
					f_(i-1,j) = u_(i-1,j);
					u_(i,j) = -u_(i,j+1);
					f_(i,j) = u_(i,j);
					break;
				case 1: //edge upper right
					u_(i,j) = 0;
					f_(i,j) = u_(i,j);
					break;
			}
		}
		else
		{
			std::cout << "Unterscheidung der Fälle für das Setzen der Obstakle-Werte u lief falsch: i: " << i << " j: " << j << std::endl;
		}
		
	}
}

void StaggeredGrid::setObstacleValues_v_g(std::array<int, 2> locations_boundary, int i, int j)
{
	// sort locations_boundary from big to small
	if(locations_boundary[0] < locations_boundary[1])
	{
		int temp = locations_boundary[0];
		locations_boundary[0] = locations_boundary[1];
		locations_boundary[1] = temp;
	}

	// set value to nan if cell is inner obstacal cell
	if(locations_boundary[0] == -1 && locations_boundary[1] == -1)
	{
		v_(i,j) = std::nan("1");
		g_(i,j) = v_(i,j);
	}

	// proof if corner or not 
	if(locations_boundary[1] == -1)
	{
		switch(locations_boundary[0]) // 0: Left, 1: Right, 2: Lower, 3: Upper
		{
			case 0:
				v_(i,j) = -v_(i-1,j);
				g_(i,j) = v_(i,j);
				break;
			case 1:
				v_(i,j) = -v_(i+1,j);
				g_(i,j) = v_(i,j);
				break;
			case 2:
				v_(i,j-1) = 0;
				g_(i,j-1) = v_(i,j-1);
				v_(i,j) = std::nan("1");
				g_(i,j) = v_(i,j);
				break;
			case 3:
				v_(i,j) = 0;
				g_(i,j) = v_(i,j);
				break;
		}
	}
	else 
	{
		if(locations_boundary[0] == 2)
		{
			switch(locations_boundary[1])
			{
				case 0: //edge lower left
					v_(i,j-1) = 0;
					g_(i,j-1) = v_(i,j-1);
					v_(i,j) = -v_(i-1,j);
					g_(i,j) = v_(i,j);
					break;
				case 1: //edge lower right
					v_(i,j-1) = 0;
					g_(i,j-1) = v_(i,j-1);
					v_(i,j) = -v_(i+1,j);
					g_(i,j) = v_(i,j);
				break;
			}
		} 
		else if(locations_boundary[0] == 3)
		{
			switch(locations_boundary[1])
			{
				case 0: //edge upper left
					v_(i,j) = 0;
					g_(i,j) = v_(i,j);
					break;
				case 1: //edge upper right
					v_(i,j) = 0;
					g_(i,j) = v_(i,j);
					break;
			}
		}
		else
		{
			std::cout << "Unterscheidung der Fälle für das Setzen der Obstakle-Werte v lief falsch: i: " << i << " j: " << j << std::endl;
		}
		
	}
}

void StaggeredGrid::setObstacleValues_p(std::array<int, 2> locations_boundary, int i, int j)
{
	// sort locations_boundary from big to small
	if(locations_boundary[0] < locations_boundary[1])
	{
		int temp = locations_boundary[0];
		locations_boundary[0] = locations_boundary[1];
		locations_boundary[1] = temp;
	}

	// set value to nan if cell is inner obstacal cell
	if(locations_boundary[0] == -1 && locations_boundary[1] == -1)
	{
		p_(i,j) = std::nan("1");
	}

	// proof if corner or not 
	if(locations_boundary[1] == -1)
	{
		int iNeighborShift;
		int jNeighborShift;
		switch(locations_boundary[0]) // 0: Left, 1: Right, 2: Lower, 3: Upper
		{
			case 0: 
				iNeighborShift = -1;
				jNeighborShift = 0;
				break;
			case 1: 
				iNeighborShift = 1; 
				jNeighborShift = 0;
				break;
			case 2:
				iNeighborShift = 0;
				jNeighborShift = -1;
				break;
			case 3:
				iNeighborShift = 0;
				jNeighborShift = 1;
		}
		p_(i,j) = p_(i+iNeighborShift, j+jNeighborShift);

	}
	else 
	{
		int iNeighborShift0;
		int jNeighborShift0;
		int iNeighborShift1;
		int jNeighborShift1;
		switch(locations_boundary[0]) // 0: Left, 1: Right, 2: Lower, 3: Upper
		{
			case 0: 
				iNeighborShift0 = -1;
				jNeighborShift0 = 0;
				break;
			case 1: 
				iNeighborShift0 = 1; 
				jNeighborShift0 = 0;
				break;
			case 2:
				iNeighborShift0 = 0;
				jNeighborShift0 = -1;
				break;
			case 3:
				iNeighborShift0 = 0;
				jNeighborShift0 = 1;
		}
		switch(locations_boundary[1]) // 0: Left, 1: Right, 2: Lower, 3: Upper
		{
			case 0: 
				iNeighborShift1 = -1;
				jNeighborShift1 = 0;
				break;
			case 1: 
				iNeighborShift1 = 1; 
				jNeighborShift1 = 0;
				break;
			case 2:
				iNeighborShift1 = 0;
				jNeighborShift1 = -1;
				break;
			case 3:
				iNeighborShift1 = 0;
				jNeighborShift1 = 1;
		}
		p_(i,j) = 1/2*(p(i+iNeighborShift0, j+jNeighborShift0) + p(i+iNeighborShift1, j+jNeighborShift1));
	}
}

void StaggeredGrid::setObstacleValues_T(std::array<int, 2> locations_boundary, int i, int j)
{
	// sort locations_boundary from big to small
	if(locations_boundary[0] < locations_boundary[1])
	{
		int temp = locations_boundary[0];
		locations_boundary[0] = locations_boundary[1];
		locations_boundary[1] = temp;
	}

	// set value to nan if cell is inner obstacal cell
	if(locations_boundary[0] == -1 && locations_boundary[1] == -1)
	{
		T_(i,j) = std::nan("1");
	}

	// proof if corner or not 
	if(locations_boundary[1] == -1)
	{
		int iNeighborShift;
		int jNeighborShift;
		switch(locations_boundary[0]) // 0: Left, 1: Right, 2: Lower, 3: Upper
		{
			case 0: 
				iNeighborShift = -1;
				jNeighborShift = 0;
				break;
			case 1: 
				iNeighborShift = 1; 
				jNeighborShift = 0;
				break;
			case 2:
				iNeighborShift = 0;
				jNeighborShift = -1;
				break;
			case 3:
				iNeighborShift = 0;
				jNeighborShift = 1;
		}
		T_(i,j) = T_(i+iNeighborShift, j+jNeighborShift);

	}
	else 
	{
		int iNeighborShift0;
		int jNeighborShift0;
		int iNeighborShift1;
		int jNeighborShift1;
		switch(locations_boundary[0]) // 0: Left, 1: Right, 2: Lower, 3: Upper
		{
			case 0: 
				iNeighborShift0 = -1;
				jNeighborShift0 = 0;
				break;
			case 1: 
				iNeighborShift0 = 1; 
				jNeighborShift0 = 0;
				break;
			case 2:
				iNeighborShift0 = 0;
				jNeighborShift0 = -1;
				break;
			case 3:
				iNeighborShift0 = 0;
				jNeighborShift0 = 1;
		}
		switch(locations_boundary[1]) // 0: Left, 1: Right, 2: Lower, 3: Upper
		{
			case 0: 
				iNeighborShift1 = -1;
				jNeighborShift1 = 0;
				break;
			case 1: 
				iNeighborShift1 = 1; 
				jNeighborShift1 = 0;
				break;
			case 2:
				iNeighborShift1 = 0;
				jNeighborShift1 = -1;
				break;
			case 3:
				iNeighborShift1 = 0;
				jNeighborShift1 = 1;
		}
		T_(i,j) = 1/2*(T_(i+iNeighborShift0, j+jNeighborShift0) + T_(i+iNeighborShift1, j+jNeighborShift1));
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
