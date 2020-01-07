#include "StaggeredGrid.h"

// StaggeredGrid::StaggeredGrid(std::array<int,2> nCells, std::array<double,2> meshWidth) :
// 	nCells_(nCells), meshWidth_(meshWidth),
// 	u_( {nCells[0]+1, nCells[1]+2},  {0*meshWidth[0],    -0.5*meshWidth[1]}, meshWidth ),
// 	v_( {nCells[0]+2, nCells[1]+1},  {-0.5*meshWidth[0], 0*meshWidth[1]}, meshWidth ),
// 	p_( {nCells[0]+2, nCells[1]+2},  {-0.5*meshWidth[0], -0.5*meshWidth[1]}, meshWidth),
// 	f_( {nCells[0]+1, nCells[1]+2},  {0*meshWidth[0],    -0.5*meshWidth[1]}, meshWidth ),
// 	g_( {nCells[0]+2, nCells[1]+1},  {-0.5*meshWidth[0], 0*meshWidth[1]}, meshWidth ),
// 	rhs_( {nCells[0]+2, nCells[1]+2},{-0.5*meshWidth[0], -0.5*meshWidth[1]}, meshWidth),
// 	T_( {nCells[0]+2, nCells[1]+2},  {-0.5*meshWidth[0], -0.5*meshWidth[1]}, meshWidth)
// {};

// StaggeredGrid::StaggeredGrid(std::array< int, 2 > nCells, std::array< double, 2 > meshWidth, std::shared_ptr<Array2D> geometryPVString, std::shared_ptr<Array2D> geometryPV1, std::shared_ptr<Array2D> geometryPV2, std::shared_ptr<Array2D> geometryTString, std::shared_ptr<Array2D> geometryT1) :
// 	nCells_(nCells), meshWidth_(meshWidth),
// 	u_( {nCells[0]+1, nCells[1]+2},  {0*meshWidth[0],    -0.5*meshWidth[1]}, meshWidth ),
// 	v_( {nCells[0]+2, nCells[1]+1},  {-0.5*meshWidth[0], 0*meshWidth[1]}, meshWidth ),
// 	p_( {nCells[0]+2, nCells[1]+2},  {-0.5*meshWidth[0], -0.5*meshWidth[1]}, meshWidth),
// 	f_( {nCells[0]+1, nCells[1]+2},  {0*meshWidth[0],    -0.5*meshWidth[1]}, meshWidth ),
// 	g_( {nCells[0]+2, nCells[1]+1},  {-0.5*meshWidth[0], 0*meshWidth[1]}, meshWidth ),
// 	rhs_( {nCells[0]+2, nCells[1]+2},{-0.5*meshWidth[0], -0.5*meshWidth[1]}, meshWidth),
// 	T_( {nCells[0]+2, nCells[1]+2},  {-0.5*meshWidth[0], -0.5*meshWidth[1]}, meshWidth),
// 	geometryPVString_(geometryPVString), geometryPV1_(geometryPV1), geometryPV2_(geometryPV2), geometryTString_(geometryTString), geometryT1_(geometryT1)
// {};

StaggeredGrid::StaggeredGrid(std::array< int, 2 > nCells, std::array< double, 2 > meshWidth, std::shared_ptr<Array2D> geometryPVString, std::shared_ptr<Array2D> geometryPVOrientation, std::shared_ptr<Array2D> geometryPV1, std::shared_ptr<Array2D> geometryPV2, std::shared_ptr<Array2D> geometryTString, std::shared_ptr<Array2D> geometryT1, Adapter& adapter) :
	nCells_(nCells), meshWidth_(meshWidth),
	u_( {nCells[0]+1, nCells[1]+2},  {0*meshWidth[0],    -0.5*meshWidth[1]}, meshWidth ),
	v_( {nCells[0]+2, nCells[1]+1},  {-0.5*meshWidth[0], 0*meshWidth[1]}, meshWidth ),
	p_( {nCells[0]+2, nCells[1]+2},  {-0.5*meshWidth[0], -0.5*meshWidth[1]}, meshWidth),
	f_( {nCells[0]+1, nCells[1]+2},  {0*meshWidth[0],    -0.5*meshWidth[1]}, meshWidth ),
	g_( {nCells[0]+2, nCells[1]+1},  {-0.5*meshWidth[0], 0*meshWidth[1]}, meshWidth ),
	rhs_( {nCells[0]+2, nCells[1]+2},{-0.5*meshWidth[0], -0.5*meshWidth[1]}, meshWidth),
	T_( {nCells[0]+2, nCells[1]+2},  {-0.5*meshWidth[0], -0.5*meshWidth[1]}, meshWidth),
	T_old_( {nCells[0]+2, nCells[1]+2},  {-0.5*meshWidth[0], -0.5*meshWidth[1]}, meshWidth),
	geometryPVString_(geometryPVString), geometryPVOrientation_(geometryPVOrientation), geometryPV1_(geometryPV1), geometryPV2_(geometryPV2), geometryTString_(geometryTString), geometryT1_(geometryT1),
	adapter(adapter)
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

double& StaggeredGrid::T_old(int i, int j)
{
  return StaggeredGrid::T_old_(i,j);
};

double StaggeredGrid::T_old(int i, int j) const
{
  return StaggeredGrid::T_old_(i,j);
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

void StaggeredGrid::saveOldStateT()
{
	for (int i = 0; i < T_.size().at(0); i++)
	{
		for (int j = 0; j < T_.size().at(1); j++)
		{
			T_old(i,j) = T(i,j);
		}
	}
}

void StaggeredGrid::reloadOldStateT()
{
	for (int i = 0; i < T_.size().at(0); i++)
	{
		for (int j = 0; j < T_.size().at(1); j++)
		{
			T(i,j) = T_old(i,j);
		}
	}
}

void StaggeredGrid::setBoundaryValues_u_f()
{
	for (int i = uIBegin(); i < uIEnd(); i++) // double check!
	{
		for (int j = uJBegin()-1; j < uJEnd()+1; j++)
		{
			// indices in geometry file (shifted by uIBegin and increased by one at the right (u) and upper(v) boundaries)
			int igeom = i-uIBegin()+1; 
			int jgeom = j-uJBegin()+1;

			int boundary_type = geometryPVString_->operator()(igeom, jgeom);

			if (boundary_type != -1)
			{ 
				// indices of the neigboring cell (can be fluid or solid)
				int in = i;
				int jn = j;

				int orientation = int(geometryPVOrientation_->operator()(igeom,jgeom));
				switch (orientation)
				{
					case 0: std::cout << "No orientation assigned!" << std::endl; break;
					case 1: in = i-1; break; // left is fluid
					case 2: in = i+1; break; // right
					case 3: jn = j-1; break; // lower
					case 4: jn = j+1; break; // upper
					case 5: in = i-1; jn = j-1; break; // lower-left
					case 6: in = i-1; jn = j+1; break; // upper-left
					case 7: in = i+1; jn = j-1; break; // lower-right
					case 8: in = i+1; jn = j+1; break; // upper-right
					default: std::cout << "Unkown orientation!" << std::endl; break;
				}

				switch (boundary_type)
				{
					case 0: // NOSLIP
					{
						switch (orientation)
						{
						case 1: case 2:
							u(i,j) = 0;
							break;
						case 3: case 4:
							u(i,j) = - u(in,jn);
							break;
						case 5: case  6: case 7: case 8:
							u(i,j) = -.5*u(i,jn);
							break;
						}
						f(i,j) = u(i,j);
						break;
					}
					case 1: // SLIP
					{
						switch (orientation)
						{
						case 1: case 2:
							u(i,j) = 0;
							break;
						case 3: case 4:
							u(i,j) = u(in,jn);
							break;
						case 5: case 6: case 7: case 8:
							u(i,j) = .5*u(i,jn);
							break;					
						}
						f(i,j) = u(i,j);
						break;
					}
					case 2: // INFLOW
					{
						switch (orientation)
						{
						case 1: case 2:
							u(i,j) = geometryPV1_->operator()(igeom,jgeom);
							break;
						case 3: case 4:
							u(i,j) = 2*geometryPV1_->operator()(igeom,jgeom) - u(in,jn);
							f(i,j) = u(i,j);
							break;
						case 5: case 6: case 7: case 8:
							u(i,j) = 3/2*geometryPV1_->operator()(igeom,jgeom) - .5*u(i,jn);
							f(i,j) = u(i,j);
							break;
						}
						f(i,j) = u(i,j);
						break;
					}
					case 3: // OUTFLOW (same as case 4)
					case 4: // PRESSURE
					{
						double u_old = u(i,j);
						switch (orientation)
						{
						case 1: case 2: case 3: case 4:
							u(i,j) = u(in,jn);
							break;
						case 5: case 6: case 7: case 8:
							u(i,j) = .5*( u(i,jn) + u(in,j) );
							break;
						}
						f(i,j) = 2*u(i,j)-u_old;
						break;
					}
				}
			}
		}
	}
}

void StaggeredGrid::setBoundaryValues_v_g()
{
	for (int i = vIBegin(); i < vIEnd()+1; i++) // double check!
	{
		for (int j = vJBegin()-1; j < vJEnd(); j++)
		{
			// indices in geometry file (shifted by uIBegin and increased by one at the right (u) and upper(v) boundaries)
			int igeom = i-vIBegin()+1; 
			int jgeom = j-vJBegin()+1;

			int boundary_type = geometryPVString_->operator()(igeom, jgeom);

			if (boundary_type != -1)
			{ 
				// indices of the neigboring cell (can be fluid or solid)
				int in = i;
				int jn = j;

				int orientation = int(geometryPVOrientation_->operator()(igeom,jgeom));
				switch (orientation)
				{
					case 0: std::cout << "No orientation assigned!" << std::endl; break;
					case 1: in = i-1; break; // left is fluid
					case 2: in = i+1; break; // right
					case 3: jn = j-1; break; // lower
					case 4: jn = j+1; break; // upper
					case 5: in = i-1; jn = j-1; break; // lower-left
					case 6: in = i-1; jn = j+1; break; // upper-left
					case 7: in = i+1; jn = j-1; break; // lower-right
					case 8: in = i+1; jn = j+1; break; // upper-right
					default: std::cout << "Unkown orientation!" << std::endl; break;
				}

				switch (boundary_type)
				{
					case 0: // NOSLIP
					{
						switch (orientation)
						{
						case 1: case 2:
							v(i,j) = - v(in,jn);
							break;
						case 3: case 4:
							v(i,j) = 0;
							break;
						case 5: case  6: case 7: case 8:
							v(i,j) = -.5*v(in,j);
							break;
						}
						g(i,j) = v(i,j);
						break;
					}
					case 1: // SLIP
					{
						switch (orientation)
						{
						case 1: case 2:
							v(i,j) = v(in,jn);
							break;
						case 3: case 4:
							v(i,j) = 0;
							break;
						case 5: case 6: case 7: case 8:
							v(i,j) = .5*v(in,j);
							break;					
						}
						g(i,j) = v(i,j);
						break;
					}
					case 2: // INFLOW
					{
						switch (orientation)
						{
						case 1: case 2:
							v(i,j) = 2*geometryPV2_->operator()(igeom,jgeom) - v(in,jn);
							g(i,j) = v(i,j);
							break;
						case 3: case 4:
							v(i,j) = geometryPV2_->operator()(igeom,jgeom);
							break;
						case 5: case 6: case 7: case 8:
							v(i,j) = 3/2*geometryPV2_->operator()(igeom,jgeom) - .5*v(in,j);
							g(i,j) = v(i,j);
							break;
						}
						g(i,j) = v(i,j);
						break;
					}
					case 3: // OUTFLOW (same as case 4)
					case 4: // PRESSURE
					{
						double u_old = v(i,j);
						switch (orientation)
						{
						case 1: case 2: case 3: case 4:
							v(i,j) = v(in,jn);
							break;
						case 5: case 6: case 7: case 8:
							v(i,j) = .5*( v(i,jn) + v(in,j) );
							break;
						}
						g(i,j) = 2*v(i,j)-u_old;
						break;
					}
				}
			}
		}
	}
}

void StaggeredGrid::setBoundaryValues_p()
{
	for (int i = pIBegin(); i < pIEnd()+1; i++) // double check!
	{
		for (int j = pJBegin()-1; j < pJEnd()+1; j++)
		{
			// indices in geometry file (shifted by uIBegin and increased by one at the right (u) and upper(v) boundaries)
			int igeom = i-pIBegin()+1; 
			int jgeom = j-pJBegin()+1;

			int boundary_type = geometryPVString_->operator()(igeom, jgeom);

			if (boundary_type != -1)
			{ 
				// indices of the neigboring cell (can be fluid or solid)
				int in = i;
				int jn = j;

				int orientation = int(geometryPVOrientation_->operator()(igeom,jgeom));
				switch (orientation)
				{
					case 0: std::cout << "No orientation assigned!" << std::endl; break;
					case 1: in = i-1; break; // left is fluid
					case 2: in = i+1; break; // right
					case 3: jn = j-1; break; // lower
					case 4: jn = j+1; break; // upper
					case 5: in = i-1; jn = j-1; break; // lower-left
					case 6: in = i-1; jn = j+1; break; // upper-left
					case 7: in = i+1; jn = j-1; break; // lower-right
					case 8: in = i+1; jn = j+1; break; // upper-right
					default: std::cout << "Unkown orientation!" << std::endl; break;
				}

				if (boundary_type == 4) // PRESSURE
				{
					switch (orientation)
					{
						case 1: case 2: case 3: case 4:
							p(i,j) = 2*geometryPV1_->operator()(igeom, jgeom) - p(in,jn);
							break;
						case 5: case 6: case 7: case 8:
							p(i,j) = 2*geometryPV1_->operator()(igeom, jgeom) - .5*(p(i,jn)+ p(in,j));
							break;
					}
				}
				else
				{
					switch (orientation)
					{
						case 1: case 2: case 3: case 4:
							p(i,j) = p(in,jn);
							break;
						case 5: case 6: case 7: case 8:
							p(i,j) = .5*(p(i,jn)+p(in,j));
							break;
					}
				}
				
			}
		}
	}
}

void StaggeredGrid::setBoundaryValues_T(double* readData, int vertexSize, std::vector<int> &vertex_i, std::vector<int> &vertex_j)
{
	for (int i = pIBegin(); i < pIEnd()+1; i++) // double check!
	{
		for (int j = pJBegin()-1; j < pJEnd()+1; j++)
		{
			// indices in geometry file (shifted by uIBegin and increased by one at the right (u) and upper(v) boundaries)
			int igeom = i-pIBegin()+1; 
			int jgeom = j-pJBegin()+1;

			int boundary_type = geometryPVString_->operator()(igeom, jgeom);

			if (boundary_type != -1)
			{ 
				// indices of the neigboring cell (can be fluid or solid)
				int in = i;
				int jn = j;

				int orientation = int(geometryPVOrientation_->operator()(igeom,jgeom));
				switch (orientation)
				{
					case 0: std::cout << "No orientation assigned!" << std::endl; break;
					case 1: in = i-1; break; // left is fluid
					case 2: in = i+1; break; // right
					case 3: jn = j-1; break; // lower
					case 4: jn = j+1; break; // upper
					case 5: in = i-1; jn = j-1; break; // lower-left
					case 6: in = i-1; jn = j+1; break; // upper-left
					case 7: in = i+1; jn = j-1; break; // lower-right
					case 8: in = i+1; jn = j+1; break; // upper-right
					default: std::cout << "Unkown orientation!" << std::endl; break;
				}

				// set boundary values of one cell
				boundary_type = geometryTString_->operator()(igeom, jgeom);
				switch (boundary_type)
				{
					case 0: // TD
						switch (orientation)
						{
							case 1: case 2: case 3: case 4:
								T(i,j) = 2*geometryT1_->operator()(igeom, jgeom) - T(in,jn);
								break;
							case 5: case 6: case 7: case 8: 
								T(i,j) = 2*geometryT1_->operator()(igeom, jgeom) - .5*(T(i,jn) + T(in,j));
								break;
						}
						break;
					case 1: // TN
						switch (orientation)
						{
							case 1: case 2: 
								T(i,j) = T(in,jn) + dx()*geometryT1_->operator()(igeom, jgeom);
								break;
							case 3: case 4:
								T(i,j) = T(in,jn) + dy()*geometryT1_->operator()(igeom, jgeom);
								break;
							case 5: case 6: case 7: case 8:
								T(i,j) = .5*(T(i,jn) + T(in,j)+ (dx()+dy())*geometryT1_->operator()(igeom, jgeom));
								break;
						}
					case 2: // TPD
					case 3: // TPN
						break;
					default: std::cout << "Unknow temperature boundary condition!" << std::endl;
				}
			}
		}
	}
	// TPD/TPN
	for (int v = 0; v < vertexSize; v++)
	{
		T(vertex_i.at(v),vertex_j.at(v)) = readData[v];
	}
}

// void StaggeredGrid::setBoundaryValues_u_f(int location_boundary, int i, int j)
// {
// 	// indices in geometry file (shifted by uIBegin and increased by one at the right (u) and upper(v) boundaries)
// 	int igeom = i-uIBegin()+1; 
// 	int jgeom = j-uJBegin()+1;
// 	// indices of the neigboring cell (can be fluid or solid)
// 	int in = i;
// 	int jn = j;
// 	switch (location_boundary) // 0: Left, 1: Right, 2: Lower, 3: Upper
// 	{
// 		case 0:	in = i+1; igeom = 0; break;
// 		case 1: in = i-1; igeom = geometryPVString_->size()[0]-1; break;
// 		case 2: jn = j+1; jgeom = 0; break;
// 		case 3: jn = j-1; jgeom = geometryPVString_->size()[1]-1; break;
// 	}
// 	// set boundary values to nan if not needed (neighbor not fluid cell)
// 	if (false){//((geometryPVString_->operator()(igeom + in - i, jgeom + jn - j) != -1)	{
// 		// std::cout << "nan1: " << geometryPVString_->operator()(igeom + in - i, jgeom + jn - j) << " " << igeom + in - i << " " << jgeom + jn - j << " " << igeom << i << in << " " << jgeom << j << jn << std::endl;
// 		// u(i,j) = std::nan("1"); commented out by Henrik
// 	}
// 	else
// 	{
// 		// set boundary values of one cell
// 		int boundary_type = geometryPVString_->operator()(igeom, jgeom);
// 		switch (boundary_type)
// 		{
// 			case 0: // NOSLIP
// 			{
// 				if (location_boundary == 0 || location_boundary == 1) { // left or right
// 					u(i,j) = 0;
// 					f(i,j) = u(i,j);
// 				} else { // upper or lower  skript meint einfach wir kompliziert
// 					u(i,j) = - u(in,jn);
// 					f(i,j) = u(i,j); //40000; //u(i,j);
// 				}
// 				break;
// 			}
// 			case 1: // SLIP
// 			{
// 				if (location_boundary == 0 || location_boundary == 1) { // left or right
// 					u(i,j) = 0;
// 					f(i,j) = u(i,j);
// 				} else { // upper or lower
// 					double u_old = u(i,j);
// 					u(i,j) = u(in,jn);
// 					f(i,j) = u(i,j); //40000; //2*u(i,j)-u_old;
// 				}
// 				break;
// 			}
// 			case 2: // INFLOW
// 			{
// 				if (location_boundary == 0 || location_boundary == 1) { // left or right
// 					u(i,j) = geometryPV1_->operator()(igeom,jgeom);
// 					f(i,j) = u(i,j);
// 				} else { // upper or lower
// 					u(i,j) = 2*geometryPV1_->operator()(igeom,jgeom) - u(in,jn);
// 					f(i,j) = u(i,j);
// 				}
// 				break;
// 			}
// 			case 3: // OUTFLOW (same as case 4)
// 			case 4: // PRESSURE
// 			{
// 				if (location_boundary == 0 || location_boundary == 1) { // left or right
// 					double u_old = u(i,j);
// 					u(i,j) = u(in,jn);
// 					f(i,j) = 2*u(i,j)-u_old;
// 				} else { // upper or lower
// 					double u_old = u(i,j);
// 					u(i,j) = u(in,jn);
// 					f(i,j) = 2*u(i,j)-u_old;
// 				}
// 				break;
// 			}

// 		}
// 	}
// 	// f(i,j) = u(i,j);

// }

// void StaggeredGrid::setBoundaryValues_v_g(int location_boundary, int i, int j)
// {
// 	// indices in geometry file (shifted by uIBegin and increased by one at the right (u) and upper(v) boundaries)
// 	int igeom = i-vIBegin()+1; // todo: double check!!
// 	int jgeom = j-vJBegin()+1;
// 	// indices of the neigboring cell (can be fluid or solid)
// 	int in = i;
// 	int jn = j;
// 	switch (location_boundary) // 0: Left, 1: Right, 2: Lower, 3: Upper
// 	{
// 		case 0:	in = i+1; igeom = 0; break;
// 		case 1: in = i-1; igeom = geometryPVString_->size()[0]-1; break;
// 		case 2: jn = j+1; jgeom = 0; break;
// 		case 3: jn = j-1; jgeom = geometryPVString_->size()[1]-1; break;
// 	}
// 	// set boundary values to nan if not needed (neighbor not fluid cell)
// 	if (false){//(geometryPVString_->operator()(igeom + in - i, jgeom + jn - j) != -1)	{
// 		// v(i,j) = std::nan("1"); commented out by Henrik
// 	}
// 	else
// 	{
// 		// set boundary values of one cell
// 		int boundary_type = geometryPVString_->operator()(igeom, jgeom);
// 		switch (boundary_type)
// 		{
// 			case 0: // NOSLIP
// 			{
// 				if (location_boundary == 0 || location_boundary == 1) { // left or right
// 					v(i,j) = - v(in,jn);
// 					g(i,j) = v(i,j);//40000; //v(i,j);
// 				} else { // upper or lower
// 					v(i,j) = 0;
// 					g(i,j) = v(i,j);
// 				}
// 				break;
// 			}
// 			case 1: // SLIP
// 			{
// 				if (location_boundary == 0 || location_boundary == 1) { // left or right
// 					double v_old = v(i,j);
// 					v(i,j) = v(in,jn);
// 					g(i,j) = v(i,j);//40000; //2*v(i,j)-v_old;
// 				} else { // upper or lower
// 					v(i,j) = 0;
// 					g(i,j) = v(i,j);
// 				}
// 				break;
// 			}
// 			case 2: // INFLOW
// 			{
// 				if (location_boundary == 0 || location_boundary == 1) { // left or right
// 					v(i,j) = 2*geometryPV2_->operator()(igeom,jgeom) - v(in,jn);
// 					g(i,j) = v(i,j);
// 				} else { // upper or lower
// 					v(i,j) = geometryPV2_->operator()(igeom,jgeom);
// 					g(i,j) = v(i,j);
// 				}
// 				break;
// 			}
// 			case 3: // OUTFLOW (same as case 4)
// 			case 4: // PRESSURE
// 			{
// 				if (location_boundary == 0 || location_boundary == 1) { // left or right
// 					double v_old = v(i,j);
// 					v(i,j) = v(in,jn);
// 					g(i,j) = 2*v(i,j)-v_old;
// 				} else { // upper or lower
// 					double v_old = v(i,j);
// 					v(i,j) = v(in,jn);
// 					g(i,j) = 2*v(i,j)-v_old;
// 				}
// 				break;
// 			}

// 		}
// 	}
// 	// g(i,j) = v(i,j);

// }

// void StaggeredGrid::setBoundaryValues_p(int location_boundary, int i, int j)
// {
// 	// indices in geometry file (shifted by uIBegin and increased by one at the right (u) and upper(v) boundaries)
// 	int igeom = i-pIBegin()+1; // todo: double check!!
// 	int jgeom = j-pJBegin()+1;
// 	// indices of the neigboring cell (can be fluid or solid)
// 	int in = i;
// 	int jn = j;
// 	switch (location_boundary) // 0: Left, 1: Right, 2: Lower, 3: Upper
// 	{
// 		case 0:	in = i+1; igeom = 0; break;
// 		case 1: in = i-1; igeom = geometryPVString_->size()[0]-1; break;
// 		case 2: jn = j+1; jgeom = 0; break;
// 		case 3: jn = j-1; jgeom = geometryPVString_->size()[1]-1; break;
// 	}
// 	// set boundary values to nan if not needed (neighbor not fluid cell)
// 	if (false){//((geometryPVString_->operator()(igeom + in - i, jgeom + jn - j) != -1)	{
// 		// std::cout << "nan3: " << location_boundary << ", " << in << ", " << i << ", " << igeom << ", j: " << jn << j << jgeom << std::endl;
// 		// p(i,j) = std::nan("1"); commented out by Henrik
// 	}
// 	else
// 	{
// 		// set boundary values of one cell
// 		int boundary_type = geometryPVString_->operator()(igeom, jgeom);
// 		if (boundary_type == 4) // PRESSURE
// 		{
// 			p(i,j) = 2*geometryPV1_->operator()(igeom, jgeom) - p(in,jn);
// 		}
// 		else
// 		{
// 			p(i,j) = p(in,jn);
// 		}

// 	}
// }

// void StaggeredGrid::setBoundaryValues_T(int location_boundary, int i, int j)
// {
// 	// indices in geometry file (shifted by uIBegin and increased by one at the right (u) and upper(v) boundaries)
// 	int igeom = i-pIBegin()+1; // todo: double check!!
// 	int jgeom = j-pJBegin()+1;
// 	// mesh width h
// 	double dx = meshWidth_[0];
// 	double dy = meshWidth_[1];
// 	// indices of the neigboring cell (can be fluid or solid)
// 	int in = i;
// 	int jn = j;
// 	// To do: choose neighbor cell like for the obstacle cells
// 	switch (location_boundary) // 0: Left, 1: Right, 2: Lower, 3: Upper
// 	{
// 		case 0:	in = i+1; igeom = 0; break;
// 		case 1: in = i-1; igeom = geometryPVString_->size()[0]-1; break;
// 		case 2: jn = j+1; jgeom = 0; break;
// 		case 3: jn = j-1; jgeom = geometryPVString_->size()[1]-1; break;
// 	}
// 	// set boundary values to nan if not needed (neighbor not fluid cell)
// 	if (false){//((geometryPVString_->operator()(igeom + in - i, jgeom + jn - j) != -1)	{
// 		// T(i,j) = std::nan("1"); commented out by Henrik
// 	}
// 	else
// 	{
// 		// set boundary values of one cell
// 		int boundary_type = geometryTString_->operator()(igeom, jgeom);
// 		if (boundary_type == 0) // TD
// 		{
// 			T(i,j) = 2*geometryT1_->operator()(igeom, jgeom) - T(in,jn);
// 		}
// 		else if (boundary_type == 1)  // TN
// 		{
// 			if(location_boundary <= 1)
// 			{
// 				T(i,j) = T(in,jn) + dx*geometryT1_->operator()(igeom, jgeom);
// 			}
// 			else
// 			{
// 				T(i,j) = T(in,jn) + dy*geometryT1_->operator()(igeom, jgeom);
// 			}
// 		}
// 		else if (boundary_type == 2) // TPD
// 		{
// 			// To do: ...
// 		} 
// 		else if (boundary_type == 3) // TPN
// 		{
// 			// To do: ...
// 		} 
// 		else 
// 		{
// 			std::cout << "Unknow temperature boundary condition!" << std::endl;
// 		}

// 	}
// }

void StaggeredGrid::setObstacleValues_u_f2(int i, int j)
{
	int igeom = i-uIBegin()+1; // todo: double check!!
	int jgeom = j-uJBegin()+1;

	if (geometryPVString_->operator()(igeom + 1, jgeom) == -1) // right is fluid
	{
		u_(i,j) = 0;
	}
	else
	{
		if (geometryPVString_->operator()(igeom - 1, jgeom) == -1) // left is fluid
		{
			u_(i-1,j) = 0;
			f_(i-1,j) = u_(i-1,j);
		}
		if (geometryPVString_->operator()(igeom, jgeom - 1) == -1 && geometryPVString_->operator()(igeom + 1, jgeom - 1) == -1) // lower and lowerright is fluid //todo: correct grammar
		{
			u_(i,j) = -u_(i,j-1);
		}
		else if (geometryPVString_->operator()(igeom, jgeom + 1) == -1 && geometryPVString_->operator()(igeom + 1, jgeom + 1) == -1) // upper and upperright is fluid
		{
			u_(i,j) = -u_(i,j+1);
		}
		else // only left is fluid and other cases
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

	if (geometryPVString_->operator()(igeom, jgeom + 1) == -1) // upper is fluid
	{
		v_(i,j) = 0;
	}
	else
	{
		if (geometryPVString_->operator()(igeom, jgeom - 1) == -1) // lower is fluid
		{
			v_(i,j-1) = 0;
			g_(i,j-1) = v_(i,j-1);
		}
		if (geometryPVString_->operator()(igeom - 1, jgeom) == -1 && geometryPVString_->operator()(igeom - 1, jgeom + 1) == -1) // left and leftupper is fluid //todo: correct grammar
		{
			v_(i,j) = -v_(i-1,j);
		}
		else if (geometryPVString_->operator()(igeom+ 1, jgeom) == -1 && geometryPVString_->operator()(igeom + 1, jgeom + 1) == -1) // right and rightupper is fluid
		{
			v_(i,j) = -v_(i+1,j);
		}
		else // only left is fluid and other cases
		{
			v_(i,j) = std::nan("1");
		}
	}
	g_(i,j) = v_(i,j);

}

void StaggeredGrid::setObstacleValues_p2(int i, int j)
{
	int igeom = i-pIBegin()+1; // todo: double check!!
	int jgeom = j-pJBegin()+1;

	if (geometryPVString_->operator()(igeom, jgeom) != 5) // if current cell is not solid
	{
		return;
	}

	double average = 0;
	int countFluidCells = 0;
	// left
	if(geometryPVString_->operator()(igeom-1, jgeom) == -1)
	{
		average += p_(i-1,j);
		countFluidCells ++;
	}

	// right
	if(geometryPVString_->operator()(igeom + 1, jgeom) == -1)
	{
		average += p_(i+1,j);
		countFluidCells++;
	}

	// below
	if(geometryPVString_->operator()(igeom, jgeom-1) == -1)
	{
		average += p_(i,j-1);
		countFluidCells++;
	}

	// upper
	if(geometryPVString_->operator()(igeom, jgeom +1) == -1)
	{
		average += p_(i,j+1);
		countFluidCells++;
	}

	if(countFluidCells > 0)
	{
		p_(i,j) = average/countFluidCells;
	}
	else
	{
		p_(i,j) = std::nan("1");
	}

}

void StaggeredGrid::setObstacleValues_T2(int i, int j)
{

	int igeom = i-pIBegin()+1; // todo: double check!!
	int jgeom = j-pJBegin()+1;

	if (geometryPVString_->operator()(igeom, jgeom) != 5) // if current cell is not solid
	{
		return;
	}

	double average = 0;
	int countFluidCells = 0;
	// left
	if(geometryPVString_->operator()(igeom-1, jgeom) == -1)
	{
		average += T_(i-1,j);
		countFluidCells ++;
	}

	// right
	if(geometryPVString_->operator()(igeom + 1, jgeom) == -1)
	{
		average += T_(i+1,j);
		countFluidCells++;
	}

	// below
	if(geometryPVString_->operator()(igeom, jgeom-1) == -1)
	{
		average += T_(i,j-1);
		countFluidCells++;
	}

	// upper
	if(geometryPVString_->operator()(igeom, jgeom +1) == -1)
	{
		average += T_(i,j+1);
		countFluidCells++;
	}

	if(countFluidCells > 0)
	{
		T_(i,j) = average/countFluidCells;
	}
	else
	{
		T_(i,j) = std::nan("1");
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


//PFUSCH!!!!!!!!!!

// // set obstaclesValues iff cell is solid for u and f
// void StaggeredGrid::setObstacleValues_u_f2(int i, int j)
// {
// 	std::array<int , 2> locations_boundary = {-2,-2};
// 	// indices in geometry file (shifted by uIBegin and increased by one at the right (u) and upper(v) boundaries)
// 	int igeom = i-uIBegin()+1; // todo: double check!!
// 	int jgeom = j-uJBegin()+1;
// 	if(geometryPVString_->operator()(igeom, jgeom) == 5) // proof, if cell is solid
// 	{
// 	int index = 0;
// 	for(int i = 0; i < 4; i++)
// 	{
// 		if(index > 2)
// 		{
// 			std::cout << "u: boundary locations are not calculated correctly or two-cell-criterion is not fullfilled, index: " << index << std::endl;
// 		}
// 		switch(i) // 0: Left, 1: Right, 2: Lower, 3: Upper
// 		{
// 			case 0:
// 				if(geometryPVString_->operator()(igeom-1, jgeom) == -1)
// 				{
// 					locations_boundary[index] = i;
// 					index++;
// 				}
// 				break;
// 			case 1:
// 				if(geometryPVString_->operator()(igeom+1, jgeom) == -1)
// 				{
// 					locations_boundary[index] = i;
// 					index++;
// 				}
// 				break;
// 			case 2:
// 				if(geometryPVString_->operator()(igeom, jgeom-1) == -1)
// 				{
// 					locations_boundary[index] = i;
// 					index++;
// 				}
// 				break;
// 			case 3:
// 				if(geometryPVString_->operator()(igeom, jgeom+1) == -1)
// 				{
// 					locations_boundary[index] = i;
// 					index++;
// 				}
// 				break;
// 		}
// 	}
//
// 	while(index < 2)
// 	{
// 		locations_boundary[index] = -1;
// 		index++;
// 	}
//
// 	setObstacleValues_u_f(locations_boundary, i, j);
// 	}
//
// }
//
// void StaggeredGrid::setObstacleValues_v_g2(int i, int j)
// {
// 	std::array<int , 2> locations_boundary = {-2,-2};
// 	// indices in geometry file (shifted by uIBegin and increased by one at the right (u) and upper(v) boundaries)
// 	int igeom = i-vIBegin()+1; // todo: double check!!
// 	int jgeom = j-vJBegin()+1;
// 	if(geometryPVString_->operator()(igeom, jgeom) == 5) // proof, if cell is solid
// 	{
// 	int index = 0;
// 	for(int i = 0; i < 4; i++)
// 	{
// 		if(index > 2)
// 		{
// 			std::cout << "v: boundary locations are not calculated correctly or two-cell-criterion is not fullfilled, index: " << index << std::endl;
// 		}
// 		switch(i) // 0: Left, 1: Right, 2: Lower, 3: Upper
// 		{
// 			case 0:
// 				if(geometryPVString_->operator()(igeom-1, jgeom) == -1)
// 				{
// 					locations_boundary[index] = i;
// 					index++;
// 				}
// 				break;
// 			case 1:
// 				if(geometryPVString_->operator()(igeom+1, jgeom) == -1)
// 				{
// 					locations_boundary[index] = i;
// 					index++;
// 				}
// 				break;
// 			case 2:
// 				if(geometryPVString_->operator()(igeom, jgeom-1) == -1)
// 				{
// 					locations_boundary[index] = i;
// 					index++;
// 				}
// 				break;
// 			case 3:
// 				if(geometryPVString_->operator()(igeom, jgeom+1) == -1)
// 				{
// 					locations_boundary[index] = i;
// 					index++;
// 				}
// 				break;
// 		}
// 	}
//
// 	while(index < 2)
// 	{
// 		locations_boundary[index] = -1;
// 		index++;
// 	}
//
// 	setObstacleValues_v_g(locations_boundary, i, j);
// 	}
// }
//
// void StaggeredGrid::setObstacleValues_p2(int i, int j)
// {
// 	std::array<int , 2> locations_boundary = {-2,-2};
// 	// indices in geometry file (shifted by uIBegin and increased by one at the right (u) and upper(v) boundaries)
// 	int igeom = i-pIBegin()+1; // todo: double check!!
// 	int jgeom = j-pJBegin()+1;
// 	if(geometryPVString_->operator()(igeom, jgeom) == 5) // proof, if cell is solid
// 	{
// 	int index = 0;
// 	for(int i = 0; i < 4; i++)
// 	{
// 		if(index > 2)
// 		{
// 			std::cout << "p: boundary locations are not calculated correctly or two-cell-criterion is not fullfilled, index: " << index << std::endl;
// 		}
// 		switch(i) // 0: Left, 1: Right, 2: Lower, 3: Upper
// 		{
// 			case 0:
// 				if(geometryPVString_->operator()(igeom-1, jgeom) == -1)
// 				{
// 					locations_boundary[index] = i;
// 					index++;
// 				}
// 				break;
// 			case 1:
// 				if(geometryPVString_->operator()(igeom+1, jgeom) == -1)
// 				{
// 					locations_boundary[index] = i;
// 					index++;
// 				}
// 				break;
// 			case 2:
// 				if(geometryPVString_->operator()(igeom, jgeom-1) == -1)
// 				{
// 					locations_boundary[index] = i;
// 					index++;
// 				}
// 				break;
// 			case 3:
// 				if(geometryPVString_->operator()(igeom, jgeom+1) == -1)
// 				{
// 					locations_boundary[index] = i;
// 					index++;
// 				}
// 				break;
// 		}
// 	}
//
// 	while(index < 2)
// 	{
// 		locations_boundary[index] = -1;
// 		index++;
// 	}
//
// 	setObstacleValues_p(locations_boundary, i, j);
// 	}
// }
//
// void StaggeredGrid::setObstacleValues_T2(int i, int j)
// {
// 	std::array<int , 2> locations_boundary = {-2,-2};
// 	// indices in geometry file (shifted by uIBegin and increased by one at the right (u) and upper(v) boundaries)
// 	int igeom = i-pIBegin()+1; // todo: double check!!
// 	int jgeom = j-pJBegin()+1;
// 	if(geometryPVString_->operator()(igeom, jgeom) == 5) // proof, if cell is solid
// 	{
// 	int index = 0;
// 	for(int i = 0; i < 4; i++)
// 	{
// 		if(index > 2)
// 		{
// 			std::cout << "T: boundary locations are not calculated correctly or two-cell-criterion is not fullfilled, index: " << index << std::endl;
// 		}
// 		switch(i) // 0: Left, 1: Right, 2: Lower, 3: Upper
// 		{
// 			case 0:
// 				if(geometryPVString_->operator()(igeom-1, jgeom) == -1)
// 				{
// 					locations_boundary[index] = i;
// 					index++;
// 				}
// 				break;
// 			case 1:
// 				if(geometryPVString_->operator()(igeom+1, jgeom) == -1)
// 				{
// 					locations_boundary[index] = i;
// 					index++;
// 				}
// 				break;
// 			case 2:
// 				if(geometryPVString_->operator()(igeom, jgeom-1) == -1)
// 				{
// 					locations_boundary[index] = i;
// 					index++;
// 				}
// 				break;
// 			case 3:
// 				if(geometryPVString_->operator()(igeom, jgeom+1) == -1)
// 				{
// 					locations_boundary[index] = i;
// 					index++;
// 				}
// 				break;
// 		}
// 	}
//
// 	while(index < 2)
// 	{
// 		locations_boundary[index] = -1;
// 		index++;
// 	}
//
// 	setObstacleValues_T(locations_boundary, i, j);
// 	}
// }
//
//
//
// void StaggeredGrid::setObstacleValues_u_f(std::array<int, 2> locations_boundary, int i, int j)
// {
// 	// sort locations_boundary from big to small
// 	if(locations_boundary[0] < locations_boundary[1])
// 	{
// 		int temp = locations_boundary[0];
// 		locations_boundary[0] = locations_boundary[1];
// 		locations_boundary[1] = temp;
// 	}
//
// 	// set value to nan if cell is inner obstacal cell
// 	if(locations_boundary[0] == -1 && locations_boundary[1] == -1)
// 	{
// 		u_(i,j) = std::nan("1");
// 		f_(i,j) = u_(i,j);
// 		return;
// 	}
//
// 	// proof if corner or not
// 	if(locations_boundary[1] == -1)
// 	{
// 		switch(locations_boundary[0]) // 0: Left, 1: Right, 2: Lower, 3: Upper
// 		{
// 			case 0:
// 				u_(i-1,j) = 0;
// 				f_(i-1,j) = u_(i-1,j);
// 				u_(i,j) = std::nan("1");
// 				f_(i,j) = u_(i,j);
// 				break;
// 			case 1:
// 				u_(i,j) = 0;
// 				f_(i,j) = u_(i,j);
// 				break;
// 			case 2:
// 				u_(i,j) = -u_(i,j-1);
// 				f_(i,j) = u_(i,j);
// 				break;
// 			case 3:
// 				u_(i,j) = -u_(i,j+1);
// 				f_(i,j) = u_(i,j);
// 				break;
// 		}
// 	}
// 	else
// 	{
// 		if(locations_boundary[0] == 2)
// 		{
// 			switch(locations_boundary[1])
// 			{
// 				case 0: //edge lower left
// 					u_(i-1,j) = 0;
// 					f_(i-1,j) = u_(i-1,j);
// 					u_(i,j) = -u_(i,j-1);
// 					f_(i,j) = u_(i,j);
// 					break;
// 				case 1: //edge lower right
// 					u_(i,j) = 0;
// 					f_(i,j) = u_(i,j);
// 					break;
// 			}
// 		}
// 		else if(locations_boundary[0] == 3)
// 		{
// 			switch(locations_boundary[1])
// 			{
// 				case 0: //edge upper left
// 					u_(i-1,j) = 0;
// 					f_(i-1,j) = u_(i-1,j);
// 					u_(i,j) = -u_(i,j+1);
// 					f_(i,j) = u_(i,j);
// 					break;
// 				case 1: //edge upper right
// 					u_(i,j) = 0;
// 					f_(i,j) = u_(i,j);
// 					break;
// 			}
// 		}
// 		else
// 		{
// 			std::cout << "Unterscheidung der Fälle für das Setzen der Obstakle-Werte u lief falsch: i: " << i << " j: " << j << std::endl;
// 		}
//
// 	}
// }
//
// void StaggeredGrid::setObstacleValues_v_g(std::array<int, 2> locations_boundary, int i, int j)
// {
// 	// sort locations_boundary from big to small
// 	if(locations_boundary[0] < locations_boundary[1])
// 	{
// 		int temp = locations_boundary[0];
// 		locations_boundary[0] = locations_boundary[1];
// 		locations_boundary[1] = temp;
// 	}
//
// 	// set value to nan if cell is inner obstacal cell
// 	if(locations_boundary[0] == -1 && locations_boundary[1] == -1)
// 	{
// 		v_(i,j) = std::nan("1");
// 		g_(i,j) = v_(i,j);
// 		return;
// 	}
//
// 	// proof if corner or not
// 	if(locations_boundary[1] == -1)
// 	{
// 		switch(locations_boundary[0]) // 0: Left, 1: Right, 2: Lower, 3: Upper
// 		{
// 			case 0:
// 				v_(i,j) = -v_(i-1,j);
// 				g_(i,j) = v_(i,j);
// 				break;
// 			case 1:
// 				v_(i,j) = -v_(i+1,j);
// 				g_(i,j) = v_(i,j);
// 				break;
// 			case 2:
// 				v_(i,j-1) = 0;
// 				g_(i,j-1) = v_(i,j-1);
// 				v_(i,j) = std::nan("1");
// 				g_(i,j) = v_(i,j);
// 				break;
// 			case 3:
// 				v_(i,j) = 0;
// 				g_(i,j) = v_(i,j);
// 				break;
// 		}
// 	}
// 	else
// 	{
// 		if(locations_boundary[0] == 2)
// 		{
// 			switch(locations_boundary[1])
// 			{
// 				case 0: //edge lower left
// 					v_(i,j-1) = 0;
// 					g_(i,j-1) = v_(i,j-1);
// 					v_(i,j) = -v_(i-1,j);
// 					g_(i,j) = v_(i,j);
// 					break;
// 				case 1: //edge lower right
// 					v_(i,j-1) = 0;
// 					g_(i,j-1) = v_(i,j-1);
// 					v_(i,j) = -v_(i+1,j);
// 					g_(i,j) = v_(i,j);
// 				break;
// 			}
// 		}
// 		else if(locations_boundary[0] == 3)
// 		{
// 			switch(locations_boundary[1])
// 			{
// 				case 0: //edge upper left
// 					v_(i,j) = 0;
// 					g_(i,j) = v_(i,j);
// 					break;
// 				case 1: //edge upper right
// 					v_(i,j) = 0;
// 					g_(i,j) = v_(i,j);
// 					break;
// 			}
// 		}
// 		else
// 		{
// 			std::cout << "Unterscheidung der Fälle für das Setzen der Obstakle-Werte v lief falsch: i: " << i << " j: " << j << std::endl;
// 		}
//
// 	}
// }
//
// void StaggeredGrid::setObstacleValues_p(std::array<int, 2> locations_boundary, int i, int j)
// {
// 	// sort locations_boundary from big to small
// 	if(locations_boundary[0] < locations_boundary[1])
// 	{
// 		int temp = locations_boundary[0];
// 		locations_boundary[0] = locations_boundary[1];
// 		locations_boundary[1] = temp;
// 	}
//
// 	// set value to nan if cell is inner obstacal cell
// 	if(locations_boundary[0] == -1 && locations_boundary[1] == -1)
// 	{
// 		p_(i,j) = std::nan("1");
// 		return;
// 	}
//
// 	// proof if corner or not
// 	if(locations_boundary[1] == -1)
// 	{
// 		int iNeighborShift;
// 		int jNeighborShift;
// 		switch(locations_boundary[0]) // 0: Left, 1: Right, 2: Lower, 3: Upper
// 		{
// 			case 0:
// 				iNeighborShift = -1;
// 				jNeighborShift = 0;
// 				break;
// 			case 1:
// 				iNeighborShift = 1;
// 				jNeighborShift = 0;
// 				break;
// 			case 2:
// 				iNeighborShift = 0;
// 				jNeighborShift = -1;
// 				break;
// 			case 3:
// 				iNeighborShift = 0;
// 				jNeighborShift = 1;
// 				break;
// 		}
// 		p_(i,j) = p_(i+iNeighborShift, j+jNeighborShift);
//
// 	}
// 	else
// 	{
// 		int iNeighborShift0;
// 		int jNeighborShift0;
// 		int iNeighborShift1;
// 		int jNeighborShift1;
// 		switch(locations_boundary[0]) // 0: Left, 1: Right, 2: Lower, 3: Upper
// 		{
// 			case 0:
// 				iNeighborShift0 = -1;
// 				jNeighborShift0 = 0;
// 				break;
// 			case 1:
// 				iNeighborShift0 = 1;
// 				jNeighborShift0 = 0;
// 				break;
// 			case 2:
// 				iNeighborShift0 = 0;
// 				jNeighborShift0 = -1;
// 				break;
// 			case 3:
// 				iNeighborShift0 = 0;
// 				jNeighborShift0 = 1;
// 				break;
// 		}
// 		switch(locations_boundary[1]) // 0: Left, 1: Right, 2: Lower, 3: Upper
// 		{
// 			case 0:
// 				iNeighborShift1 = -1;
// 				jNeighborShift1 = 0;
// 				break;
// 			case 1:
// 				iNeighborShift1 = 1;
// 				jNeighborShift1 = 0;
// 				break;
// 			case 2:
// 				iNeighborShift1 = 0;
// 				jNeighborShift1 = -1;
// 				break;
// 			case 3:
// 				iNeighborShift1 = 0;
// 				jNeighborShift1 = 1;
// 				break;
// 		}
// 		p_(i,j) = 1/2*(p(i+iNeighborShift0, j+jNeighborShift0) + p(i+iNeighborShift1, j+jNeighborShift1));
// 	}
// }
//
// void StaggeredGrid::setObstacleValues_T(std::array<int, 2> locations_boundary, int i, int j)
// {
// 	// sort locations_boundary from big to small
// 	if(locations_boundary[0] < locations_boundary[1])
// 	{
// 		int temp = locations_boundary[0];
// 		locations_boundary[0] = locations_boundary[1];
// 		locations_boundary[1] = temp;
// 	}
//
// 	// set value to nan if cell is inner obstacal cell
// 	if(locations_boundary[0] == -1 && locations_boundary[1] == -1)
// 	{
// 		T_(i,j) = std::nan("1");
// 		return;
// 	}
//
// 	// proof if corner or not
// 	if(locations_boundary[1] == -1)
// 	{
// 		int iNeighborShift;
// 		int jNeighborShift;
// 		switch(locations_boundary[0]) // 0: Left, 1: Right, 2: Lower, 3: Upper
// 		{
// 			case 0:
// 				iNeighborShift = -1;
// 				jNeighborShift = 0;
// 				break;
// 			case 1:
// 				iNeighborShift = 1;
// 				jNeighborShift = 0;
// 				break;
// 			case 2:
// 				iNeighborShift = 0;
// 				jNeighborShift = -1;
// 				break;
// 			case 3:
// 				iNeighborShift = 0;
// 				jNeighborShift = 1;
// 				break;
// 		}
// 		T_(i,j) = T_(i+iNeighborShift, j+jNeighborShift);
//
// 	}
// 	else
// 	{
// 		int iNeighborShift0;
// 		int jNeighborShift0;
// 		int iNeighborShift1;
// 		int jNeighborShift1;
// 		switch(locations_boundary[0]) // 0: Left, 1: Right, 2: Lower, 3: Upper
// 		{
// 			case 0:
// 				iNeighborShift0 = -1;
// 				jNeighborShift0 = 0;
// 				break;
// 			case 1:
// 				iNeighborShift0 = 1;
// 				jNeighborShift0 = 0;
// 				break;
// 			case 2:
// 				iNeighborShift0 = 0;
// 				jNeighborShift0 = -1;
// 				break;
// 			case 3:
// 				iNeighborShift0 = 0;
// 				jNeighborShift0 = 1;
// 				break;
// 		}
// 		switch(locations_boundary[1]) // 0: Left, 1: Right, 2: Lower, 3: Upper
// 		{
// 			case 0:
// 				iNeighborShift1 = -1;
// 				jNeighborShift1 = 0;
// 				break;
// 			case 1:
// 				iNeighborShift1 = 1;
// 				jNeighborShift1 = 0;
// 				break;
// 			case 2:
// 				iNeighborShift1 = 0;
// 				jNeighborShift1 = -1;
// 				break;
// 			case 3:
// 				iNeighborShift1 = 0;
// 				jNeighborShift1 = 1;
// 				break;
// 		}
// 		T_(i,j) = 1/2*(T_(i+iNeighborShift0, j+jNeighborShift0) + T_(i+iNeighborShift1, j+jNeighborShift1));
// 	}
//
// }
