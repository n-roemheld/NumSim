#pragma once

#include <array>
#include "FieldVariable.h"
#include "Partitioning.h"
#include <iostream>
#include <mpi.h>

class StaggeredGrid
{
public:
  //!constructor
  StaggeredGrid(std::array< int, 2 > nCells, std::array< double, 2 > meshWidth);
  StaggeredGrid(std::array< int, 2 > nCells, std::array< double, 2 > meshWidth, Partitioning parti);


  const std::array< double, 2 > 	meshWidth() const;

  const std::array< int, 2 > 	nCells() const;

  const FieldVariable & 	u() const;

  const FieldVariable & 	v() const;

  const FieldVariable & 	p() const;

  double 	u(int i, int j) const;

  double& 	u(int i, int j);

  double 	v(int i, int j) const;

  double& 	v(int i, int j);

  double 	p(int i, int j) const;

  double& 	p(int i, int j);

  double& 	rhs(int i, int j);

  double& 	f(int i, int j);

  double& 	g(int i, int j);

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

  // void set_partitioning(int MPI_rank, std::array<int,4> ranks_neighbors, std::array<bool,4> is_boundary, std::array<int,2> nCells, std::array<int,2> nCellsGlobal);

  int ownRankNo();

  bool is_boundary(int direction) const;// 0 = lower, 1 = right, 2 = upper, 3 = left

  // bool is_boundary(int direction) const;// 0 = lower, 1 = right, 2 = upper, 3 = left


  int rank_neighbor(int direction);

  std::array<int,2> nCells();

  std::array<int,2> nCellsGlobal();

  // Partitioning& partitioning();

  // const Partitioning partitioning();

  // void receive_boundary(int sender_tag, int index_fixed, int index_begin, int index_end, int source_rank, bool do_nothing, char function, bool horizontal_communication);
  // void send_boundary(int receiver_tag, int index_fixed, int index_begin, int index_end, int target_rank, bool do_nothing, char function, bool horizontal_communication);

  void              write_recv_buffer(int index_fixed, int index_begin, int index_end, char function, bool horizontal_communication, std::vector<double> rcv_buffer);
  std::vector<double> get_send_buffer(int index_fixed, int index_begin, int index_end, char function, bool horizontal_communication);

  void velocity_horizontal_communication(char u_or_f, char v_or_g);
  void   velocity_vertical_communication(char u_or_f, char v_or_g);
  void pressure_communication();


protected:
	const std::array< int, 2 > nCells_;
	const std::array< double, 2 > meshWidth_;
	FieldVariable 	u_;
	FieldVariable 	v_;
	FieldVariable 	p_;
	FieldVariable 	rhs_;
	FieldVariable 	f_;
	FieldVariable 	g_;

  Partitioning partitioning_;


};
