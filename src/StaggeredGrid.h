#pragma once

#include <array>
#include "FieldVariable.h"
#include <memory>
#include <cmath>
#include <iostream>
#include "Adapter.h"

class StaggeredGrid
{
public:
  //!constructor
  // StaggeredGrid(std::array< int, 2 > nCells, std::array< double, 2 > meshWidth);
  // StaggeredGrid(std::array< int, 2 > nCells, std::array< double, 2 > meshWidth, std::shared_ptr<Array2D> geometryPVString, std::shared_ptr<Array2D> geometryPV1, std::shared_ptr<Array2D> geometryPV2, std::shared_ptr<Array2D> geometryTString, std::shared_ptr<Array2D> geometryT1);
  StaggeredGrid(std::array< int, 2 > nCells, std::array< double, 2 > meshWidth, std::shared_ptr<Array2D> geometryPVString, std::shared_ptr<Array2D> geometryPVOrientation, std::shared_ptr<Array2D> geometryPV1, std::shared_ptr<Array2D> geometryPV2, std::shared_ptr<Array2D> geometryTString, std::shared_ptr<Array2D> geometryT1, Adapter& adapter);

  Adapter& adapter;

  const std::array< double, 2 > 	meshWidth() const;

  const std::array< int, 2 > 	nCells() const;

  const FieldVariable & 	u() const;

  const FieldVariable & 	v() const;

  const FieldVariable & 	p() const;

  const FieldVariable &   T() const;


  double 	u(int i, int j) const;

  double& 	u(int i, int j);

  double 	v(int i, int j) const;

  double& 	v(int i, int j);

  double 	p(int i, int j) const;

  double& 	p(int i, int j);

  double& 	rhs(int i, int j);

  double& 	f(int i, int j);

  double& 	g(int i, int j);

  double   T(int i, int j) const;

  double&   T(int i, int j);


  double geometryPVString(int i, int j) const;


  double 	dx() const;

  double 	dy() const;

  int 	uIBegin() const;

  int 	uIEnd() const;

  int 	uJBegin() const;

  int 	uJEnd() const;

  int 	vIBegin() const;

  int 	vIEnd() const;

  int 	vJBegin() const;

  int 	vJEnd() const;

  int 	pIBegin() const;

  int 	pIEnd() const;

  int 	pJBegin() const;

  int 	pJEnd() const;

  void saveOldStateT();
  void reloadOldStateT();

	// void setBoundaryValues_u_f(int location_boundary, int i, int j);
	// void setBoundaryValues_v_g(int location_boundary, int i, int j);
	// void setBoundaryValues_p(int location_boundary, int i, int j);
	// void setBoundaryValues_T(int location_boundary, int i, int j);
	void setBoundaryValues_u_f();
	void setBoundaryValues_v_g();
	void setBoundaryValues_p();
	void setBoundaryValues_T(double* readData, int vertexSize, 	std::vector<int> &vertex_i, 	std::vector<int> &vertex_j);

  // void setObstacleValues_u_f(int i, int j);
  // void setObstacleValues_v_g(int i, int j);
  // void setObstacleValues_p(int i, int j);
  // void setObstacleValues_T(int i, int j);

  void setObstacleValues_u_f2(int i, int j);
  void setObstacleValues_v_g2(int i, int j);
  void setObstacleValues_p2(int i, int j);
  void setObstacleValues_T2(int i, int j);

  void setObstacleValues_u_f(std::array<int, 2> locations_boundary, int i, int j);
  void setObstacleValues_v_g(std::array<int, 2> locations_boundary, int i, int j);
  void setObstacleValues_p(std::array<int, 2> locations_boundary, int i, int j);
  void setObstacleValues_T(std::array<int, 2> locations_boundary, int i, int j);

  void fillIn(int uInit, int vInit, int pInit, int TInit);

  // void setObstacleValues(int i, int j);


protected:

  double   T_old(int i, int j) const;

  double&   T_old(int i, int j);

  
	const std::array< int, 2 > nCells_;
	const std::array< double, 2 > meshWidth_;
	FieldVariable 	u_;
	FieldVariable 	v_;
	FieldVariable 	p_;
	FieldVariable 	rhs_;
	FieldVariable 	f_;
	FieldVariable 	g_;
  FieldVariable   T_;
  FieldVariable   T_old_;

  std::shared_ptr<Array2D> geometryPVString_; //< describes typ of cell for pressure
  std::shared_ptr<Array2D> geometryPVOrientation_; //< describes typ of cell for pressure
  std::shared_ptr<Array2D> geometryPV1_;
  std::shared_ptr<Array2D> geometryPV2_;

  std::shared_ptr<Array2D> geometryTString_;
  std::shared_ptr<Array2D> geometryT1_;
};
