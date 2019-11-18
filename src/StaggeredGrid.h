#pragma once

#include <array>
#include "FieldVariable.h"
#include "Partitioning.h"

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

  void send_boundary_vertical_u(int direction, int j_fixed, int i_begin, int i_end, int target_rank, bool do_nothing);
  void send_boundary_horizontal_u(int direction, int i_fixed, int j_begin, int j_end, int target_rank, bool do_nothing);
  void receive_boundary_vertical_u(MPI_Request& current_request, int sender_tag, int j_fixed, int i_begin, int i_end, int source_rank, bool do_nothing);
  void receive_boundary_horizontal_u(MPI_Request& current_request, int sender_tag, int i_fixed, int j_begin, int j_end, int source_rank, bool do_nothing);

  void send_boundary_vertical_v(int direction, int j_fixed, int i_begin, int i_end, int target_rank, bool do_nothing);
  void send_boundary_horizontal_v(int direction, int i_fixed, int j_begin, int j_end, int target_rank, bool do_nothing);
  void receive_boundary_vertical_v(MPI_Request& current_request, int sender_tag, int j_fixed, int i_begin, int i_end, int source_rank, bool do_nothing);
  void receive_boundary_horizontal_v(MPI_Request& current_request, int sender_tag, int i_fixed, int j_begin, int j_end, int source_rank, bool do_nothing);
  
  void send_boundary_vertical_f(int direction, int j_fixed, int i_begin, int i_end, int target_rank, bool do_nothing);
  void send_boundary_horizontal_f(int direction, int i_fixed, int j_begin, int j_end, int target_rank, bool do_nothing);
  void receive_boundary_vertical_f(MPI_Request& current_request, int sender_tag, int j_fixed, int i_begin, int i_end, int source_rank, bool do_nothing);
  void receive_boundary_horizontal_f(MPI_Request& current_request, int sender_tag, int i_fixed, int j_begin, int j_end, int source_rank, bool do_nothing);

  void send_boundary_vertical_g(int direction, int j_fixed, int i_begin, int i_end, int target_rank, bool do_nothing);
  void send_boundary_horizontal_g(int direction, int i_fixed, int j_begin, int j_end, int target_rank, bool do_nothing);
  void receive_boundary_vertical_g(MPI_Request& current_request, int sender_tag, int j_fixed, int i_begin, int i_end, int source_rank, bool do_nothing);
  void receive_boundary_horizontal_g(MPI_Request& current_request, int sender_tag, int i_fixed, int j_begin, int j_end, int source_rank, bool do_nothing);

  void send_boundary_vertical_p(int direction, int j_fixed, int i_begin, int i_end, int target_rank, bool do_nothing);
  void send_boundary_horizontal_p(int direction, int i_fixed, int j_begin, int j_end, int target_rank, bool do_nothing);
  void receive_boundary_vertical_p(MPI_Request& current_request, int sender_tag, int j_fixed, int i_begin, int i_end, int source_rank, bool do_nothing);
  void receive_boundary_horizontal_p(MPI_Request& current_request, int sender_tag, int i_fixed, int j_begin, int j_end, int source_rank, bool do_nothing);


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
