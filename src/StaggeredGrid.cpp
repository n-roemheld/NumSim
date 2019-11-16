#include "StaggeredGrid.h"

StaggeredGrid::StaggeredGrid(std::array<int,2> nCells, std::array<double,2> meshWidth) :
	nCells_(nCells), meshWidth_(meshWidth),
	u_( {nCells[0]+3, nCells[1]+2},  {-1*meshWidth[0],    -0.5*meshWidth[1]}, meshWidth ),
	v_( {nCells[0]+2, nCells[1]+3},  {-0.5*meshWidth[0], -1*meshWidth[1]}, meshWidth ),
	p_( {nCells[0]+2, nCells[1]+2},  {-0.5*meshWidth[0], -0.5*meshWidth[1]}, meshWidth),
	f_( {nCells[0]+3, nCells[1]+2},  {-1*meshWidth[0],    -0.5*meshWidth[1]}, meshWidth ),
	g_( {nCells[0]+2, nCells[1]+3},  {-0.5*meshWidth[0], -1*meshWidth[1]}, meshWidth ),
	rhs_( {nCells[0]+2, nCells[1]+2},{-0.5*meshWidth[0], -0.5*meshWidth[1]}, meshWidth)
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
    if (is_boundary(3)) {
        return 2; 
    }
	return 1; // one additional column to the left for inner boundaries
};

int StaggeredGrid::uIEnd() const
{
	return nCells_[0]+2; // reduced by one to account for boundary values // or not? // or not not?
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
    if (is_boundary(0)) {
        return 2;
    }
	return 1;
};

int StaggeredGrid::vJEnd() const
{
	return nCells_[1]+2; // reduced by one to account for boundary values // or not? // or not not?
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

void StaggeredGrid::set_partitioning(MPI_rank, ranks_neighbors, is_boundary, nCells, nCellsGlobal)
  {
    partitioning_ = new Partitioning(MPI_rank, ranks_neighbors, is_boundary, nCells, nCellsGlobal);
  };

int StaggeredGrid::ownRankNo()
{
    return partitioning_.ownRankNo();
};

int StaggeredGrid::rank_neighbors(int direction)
{
    return partitioning_.ranks_neighbors(direction);
};

bool StaggeredGrid::is_boundary(int direction)
{
    return partitioning_.is_boundary(direction);
};

std::array<int,2> StaggeredGrid::nCells()
{
    return partitioning_.nCells();
};

std::array<int,2> StaggeredGrid::nCellsGlobal()
{
    return partitioning_.nCellsglobal();
};

const Partitioning& partitioning()
  {
	return partitioning_;
  };

void StaggeredGrid::send_boundary_vertical(int tag, int j_fixed, int i_begin, int i_end, int target_rank, double (*fVar)(int, int), bool do_nothing)
{
    if (!do_nothing)
    {
        std::vector<double> send_buffer(i_end - i_begin);
        int bi = 0;
        for(int i = i_begin, i < i_end, i++)
        {
            send_buffer[bi] = fVar(i,j);
            bi++;
        }
        MPI_Isend(&send_buffer, i_end-i_begin, MPI_DOUBLE, target_rank, direction, MPI_COMM_WORLD);

    }
};

void StaggeredGrid::send_boundary_horizontal(int tag, int i_fixed, int j_begin, int j_end, int target_rank, double (*fVar)(int, int), bool do_nothing)
{
    if (!do_nothing)
    {
        std::vector<double> send_buffer(j_end - j_begin);
        int bj = 0;
        for(int j = j_begin, j < j_end, j++)
        {
            send_buffer[bj] = fVar(i,j);
            bj++;
        }
        MPI_Isend(&send_buffer, j_end-j_begin, MPI_DOUBLE, target_rank, direction, MPI_COMM_WORLD);

    }
};

MPI_Request StaggeredGrid::receive_boundary_vertical(int sender_tag, int j_fixed, int i_begin, int i_end, int source_rank, double (*fVar)(int, int), bool do_nothing)
{
    if (!do_nothing)
    {
        MPI_Request current_request;
        std::vector<double> rcv_buffer;
        MPI_Irecv(&rcv_buffer, i_end - i_begin, MPI_DOUBLE, source_rank, sender_tag, MPI_COMM_WORLD, &current_request);
        int j = j_fixed;
        int bi = 0;
        for (int i = i_begin; i < i_end; i++)
        {
            fVar(i,j) = rcv_buffer(bi);
            bi++;
        }
    }
}

MPI_Request StaggeredGrid::receive_boundary_horizontal(int sender_tag, int i_fixed, int j_begin, int j_end, int source_rank, double (*fVar)(int, int), bool do_nothing)
{
    if (!do_nothing)
    {
        MPI_Request current_request;
        std::vector<double> rcv_buffer;
        MPI_Irecv(&rcv_buffer, j_end - j_begin, MPI_DOUBLE, source_rank, sender_tag, MPI_COMM_WORLD, &current_request);
        int i = i_fixed;
        int bj = 0;
        for (int j = j_begin; j < j_end; j++)
        {
            fVar(i,j) = rcv_buffer(bj);
            bj++;
        }
    }
}

// old non parallel version !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// #include "StaggeredGrid.h"
//
// StaggeredGrid::StaggeredGrid(std::array<int,2> nCells, std::array<double,2> meshWidth) :
// 	nCells_(nCells), meshWidth_(meshWidth),
// 	u_( {nCells[0]+1, nCells[1]+2},  {0*meshWidth[0],    -0.5*meshWidth[1]}, meshWidth ),
// 	v_( {nCells[0]+2, nCells[1]+1},  {-0.5*meshWidth[0], 0*meshWidth[1]}, meshWidth ),
// 	p_( {nCells[0]+2, nCells[1]+2},  {-0.5*meshWidth[0], -0.5*meshWidth[1]}, meshWidth),
// 	f_( {nCells[0]+1, nCells[1]+2},  {0*meshWidth[0],    -0.5*meshWidth[1]}, meshWidth ),
// 	g_( {nCells[0]+2, nCells[1]+1},  {-0.5*meshWidth[0], 0*meshWidth[1]}, meshWidth ),
// 	rhs_( {nCells[0]+2, nCells[1]+2},{-0.5*meshWidth[0], -0.5*meshWidth[1]}, meshWidth)
// {};
//
// const std::array< double, 2 > StaggeredGrid::meshWidth() const
// {
//   return meshWidth_;
// };
//
// const std::array< int, 2 > 	StaggeredGrid::nCells() const
// {
//   return nCells_;
// };
//
// // removed const after function
// const FieldVariable& StaggeredGrid::u() const
// {
//   return u_;
// };
//
// // removed const after function
// const FieldVariable& StaggeredGrid::v() const
// {
//   return v_;
// };
//
// // removed const after function
// const FieldVariable& StaggeredGrid::p() const
// {
//   return p_;
// };
//
// double StaggeredGrid::u(int i, int j) const
// {
//   return StaggeredGrid::u_(i,j);
// };
//
// double StaggeredGrid::v(int i, int j) const
// {
//   return StaggeredGrid::v_(i,j);
// };
//
// double StaggeredGrid::p(int i, int j) const
// {
//   return StaggeredGrid::p_(i,j);
// };
//
// double& StaggeredGrid::u(int i, int j)
// {
//   return StaggeredGrid::u_(i,j);
// };
//
// double& StaggeredGrid::v(int i, int j)
// {
//   return StaggeredGrid::v_(i,j);
// };
//
// double& StaggeredGrid::p(int i, int j)
// {
//   return StaggeredGrid::p_(i,j);
// };
//
// double& StaggeredGrid::rhs(int i, int j)
// {
//   return StaggeredGrid::rhs_(i,j);
// };
//
// double& StaggeredGrid::f(int i, int j)
// {
//   return StaggeredGrid::f_(i,j);
// };
//
// double& StaggeredGrid::g(int i, int j)
// {
//   return StaggeredGrid::g_(i,j);
// };
//
// double StaggeredGrid::dx() const
// {
//   return meshWidth_[0];
// };
//
// double StaggeredGrid::dy() const
// {
//   return meshWidth_[1];
// };
//
// int StaggeredGrid::uIBegin() const
// {
// 	return 1;
// };
//
// int StaggeredGrid::uIEnd() const
// {
// 	return nCells_[0]+1; // reduced by one to account for boundary values // or not? // or not not?
// };
//
// int StaggeredGrid::uJBegin() const
// {
// 	return 1;
// };
//
// int StaggeredGrid::uJEnd() const
// {
// 	return nCells_[1]+1;
// };
//
// int StaggeredGrid::vIBegin() const
// {
// 	return 1;
// };
//
// int StaggeredGrid::vIEnd() const
// {
// 	return nCells_[0]+1;
// };
//
// int StaggeredGrid::vJBegin() const
// {
// 	return 1;
// };
//
// int StaggeredGrid::vJEnd() const
// {
// 	return nCells_[1]+1; // reduced by one to account for boundary values // or not? // or not not?
// };
//
// int StaggeredGrid::pIBegin() const
// {
// 	return 1;
// };
//
// int StaggeredGrid::pIEnd() const
// {
// 	return nCells_[0] +1;
// };
//
// int StaggeredGrid::pJBegin() const
// {
// 	return 1;
// };
//
// int StaggeredGrid::pJEnd() const
// {
// 	return nCells_[1] + 1;
// };
