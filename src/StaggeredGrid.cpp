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

StaggeredGrid::StaggeredGrid(std::array<int,2> nCells, std::array<double,2> meshWidth, Partitioning parti) :
	nCells_(nCells), meshWidth_(meshWidth),
	u_( {nCells[0]+3, nCells[1]+2},  {-1*meshWidth[0],    -0.5*meshWidth[1]}, meshWidth ),
	v_( {nCells[0]+2, nCells[1]+3},  {-0.5*meshWidth[0], -1*meshWidth[1]}, meshWidth ),
	p_( {nCells[0]+2, nCells[1]+2},  {-0.5*meshWidth[0], -0.5*meshWidth[1]}, meshWidth),
	f_( {nCells[0]+3, nCells[1]+2},  {-1*meshWidth[0],    -0.5*meshWidth[1]}, meshWidth ),
	g_( {nCells[0]+2, nCells[1]+3},  {-0.5*meshWidth[0], -1*meshWidth[1]}, meshWidth ),
	rhs_( {nCells[0]+2, nCells[1]+2},{-0.5*meshWidth[0], -0.5*meshWidth[1]}, meshWidth),
	partitioning_(parti)
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

// void StaggeredGrid::set_partitioning(int MPI_rank, std::array<int,4> ranks_neighbors, std::array<bool,4> is_boundary, std::array<int,2> nCells, std::array<int,2> nCellsGlobal)
//   {
//     partitioning_ = Partitioning(MPI_rank, ranks_neighbors, is_boundary, nCells, nCellsGlobal);
//   };

int StaggeredGrid::ownRankNo()
{
    return partitioning_.ownRankNo();
};

int StaggeredGrid::rank_neighbor(int direction)
{
    return partitioning_.rank_neighbor(direction);
};

bool StaggeredGrid::is_boundary(int direction) const
{
    return partitioning_.is_boundary(direction);
};

std::array<int,2> StaggeredGrid::nCells()
{
    return partitioning_.nCells();
};

std::array<int,2> StaggeredGrid::nCellsGlobal()
{
    return partitioning_.nCellsGlobal();
};

// Partitioning& partitioning()
//   {
// 	return partitioning_;
//   };

// const	Partitioning partitioning()
//   {
// 	return partitioning_;
//   };

void StaggeredGrid::send_boundary_vertical_u(int tag, int j_fixed, int i_begin, int i_end, int target_rank, bool do_nothing)
{
    if (!do_nothing)
    {
        std::vector<double> send_buffer(i_end - i_begin);
        int bi = 0;
        for(int i = i_begin; i < i_end; i++)
        {
            send_buffer[bi] = u(i,j_fixed);
            bi++;
        }
				MPI_Request request;
        MPI_Isend(&send_buffer, i_end-i_begin, MPI_DOUBLE, target_rank, tag, MPI_COMM_WORLD, &request);
				MPI_Request_free(&request);
    }
};
void StaggeredGrid::send_boundary_vertical_v(int tag, int j_fixed, int i_begin, int i_end, int target_rank, bool do_nothing)
{
    if (!do_nothing)
    {
        std::vector<double> send_buffer(i_end - i_begin);
        int bi = 0;
        for(int i = i_begin; i < i_end; i++)
        {
            send_buffer[bi] = v(i,j_fixed);
            bi++;
        }
				MPI_Request request;
        MPI_Isend(&send_buffer, i_end-i_begin, MPI_DOUBLE, target_rank, tag, MPI_COMM_WORLD, &request);
				MPI_Request_free(&request);
    }
};
void StaggeredGrid::send_boundary_vertical_f(int tag, int j_fixed, int i_begin, int i_end, int target_rank, bool do_nothing)
{
    if (!do_nothing)
    {
        std::vector<double> send_buffer(i_end - i_begin);
        int bi = 0;
        for(int i = i_begin; i < i_end; i++)
        {
            send_buffer[bi] = f(i,j_fixed);
            bi++;
        }
				MPI_Request request;
        MPI_Isend(&send_buffer, i_end-i_begin, MPI_DOUBLE, target_rank, tag, MPI_COMM_WORLD, &request);
				MPI_Request_free(&request);
    }
};
void StaggeredGrid::send_boundary_vertical_g(int tag, int j_fixed, int i_begin, int i_end, int target_rank, bool do_nothing)
{
    if (!do_nothing)
    {
        std::vector<double> send_buffer(i_end - i_begin);
        int bi = 0;
        for(int i = i_begin; i < i_end; i++)
        {
            send_buffer[bi] = g(i,j_fixed);
            bi++;
        }
				MPI_Request request;
        MPI_Isend(&send_buffer, i_end-i_begin, MPI_DOUBLE, target_rank, tag, MPI_COMM_WORLD, &request);
				MPI_Request_free(&request);
    }
};
void StaggeredGrid::send_boundary_vertical_p(int tag, int j_fixed, int i_begin, int i_end, int target_rank, bool do_nothing)
{
    if (!do_nothing)
    {
        std::vector<double> send_buffer(i_end - i_begin);
        int bi = 0;
        for(int i = i_begin; i < i_end; i++)
        {
            send_buffer[bi] = p(i,j_fixed);
            bi++;
        }
				MPI_Request request;
        MPI_Isend(&send_buffer, i_end-i_begin, MPI_DOUBLE, target_rank, tag, MPI_COMM_WORLD, &request);
				MPI_Request_free(&request);
    }
};

void StaggeredGrid::send_boundary_horizontal_u(int tag, int i_fixed, int j_begin, int j_end, int target_rank, bool do_nothing)
{
    if (!do_nothing)
    {
        std::vector<double> send_buffer(j_end - j_begin);
        int bj = 0;
        for(int j = j_begin; j < j_end; j++)
        {
            send_buffer[bj] = u(i_fixed,j);
            bj++;
        }
				MPI_Request request;
        MPI_Isend(&send_buffer, j_end-j_begin, MPI_DOUBLE, target_rank, tag, MPI_COMM_WORLD, &request);
				MPI_Request_free(&request);
    }
};
void StaggeredGrid::send_boundary_horizontal_v(int tag, int i_fixed, int j_begin, int j_end, int target_rank, bool do_nothing)
{
    if (!do_nothing)
    {
        std::vector<double> send_buffer(j_end - j_begin);
        int bj = 0;
        for(int j = j_begin; j < j_end; j++)
        {
            send_buffer[bj] = v(i_fixed,j);
            bj++;
        }
				MPI_Request request;
        MPI_Isend(&send_buffer, j_end-j_begin, MPI_DOUBLE, target_rank, tag, MPI_COMM_WORLD, &request);
				MPI_Request_free(&request);
    }
};void StaggeredGrid::send_boundary_horizontal_f(int tag, int i_fixed, int j_begin, int j_end, int target_rank, bool do_nothing)
{
    if (!do_nothing)
    {
        std::vector<double> send_buffer(j_end - j_begin);
        int bj = 0;
        for(int j = j_begin; j < j_end; j++)
        {
            send_buffer[bj] = f(i_fixed,j);
            bj++;
        }
				MPI_Request request;
        MPI_Isend(&send_buffer, j_end-j_begin, MPI_DOUBLE, target_rank, tag, MPI_COMM_WORLD, &request);
				MPI_Request_free(&request);
    }
};void StaggeredGrid::send_boundary_horizontal_g(int tag, int i_fixed, int j_begin, int j_end, int target_rank, bool do_nothing)
{
    if (!do_nothing)
    {
        std::vector<double> send_buffer(j_end - j_begin);
        int bj = 0;
        for(int j = j_begin; j < j_end; j++)
        {
            send_buffer[bj] = g(i_fixed,j);
            bj++;
        }
				MPI_Request request;
        MPI_Isend(&send_buffer, j_end-j_begin, MPI_DOUBLE, target_rank, tag, MPI_COMM_WORLD, &request);
				MPI_Request_free(&request);
    }
};void StaggeredGrid::send_boundary_horizontal_p(int tag, int i_fixed, int j_begin, int j_end, int target_rank, bool do_nothing)
{
    if (!do_nothing)
    {
        std::vector<double> send_buffer(j_end - j_begin);
        int bj = 0;
        for(int j = j_begin; j < j_end; j++)
        {
            send_buffer[bj] = g(i_fixed,j);
            bj++;
        }
				MPI_Request request;
        MPI_Isend(&send_buffer, j_end-j_begin, MPI_DOUBLE, target_rank, tag, MPI_COMM_WORLD, &request);
				MPI_Request_free(&request);
    }
};


MPI_Request StaggeredGrid::receive_boundary_vertical_u(int sender_tag, int j_fixed, int i_begin, int i_end, int source_rank, bool do_nothing)
{

		MPI_Request current_request;
    if (!do_nothing)
    {
        std::vector<double> rcv_buffer;
        MPI_Irecv(&rcv_buffer, i_end - i_begin, MPI_DOUBLE, source_rank, sender_tag, MPI_COMM_WORLD, &current_request);
        int bi = 0;
        for (int i = i_begin; i < i_end; i++)
        {
            u(i,j_fixed) = rcv_buffer[bi];
            bi++;
        }
    }
		return current_request;
}
MPI_Request StaggeredGrid::receive_boundary_vertical_v(int sender_tag, int j_fixed, int i_begin, int i_end, int source_rank, bool do_nothing)
{

		MPI_Request current_request;
    if (!do_nothing)
    {
        std::vector<double> rcv_buffer;
        MPI_Irecv(&rcv_buffer, i_end - i_begin, MPI_DOUBLE, source_rank, sender_tag, MPI_COMM_WORLD, &current_request);
        int bi = 0;
        for (int i = i_begin; i < i_end; i++)
        {
            v(i,j_fixed) = rcv_buffer[bi];
            bi++;
        }
    }
		return current_request;
}
MPI_Request StaggeredGrid::receive_boundary_vertical_f(int sender_tag, int j_fixed, int i_begin, int i_end, int source_rank, bool do_nothing)
{

		MPI_Request current_request;
    if (!do_nothing)
    {
        std::vector<double> rcv_buffer;
        MPI_Irecv(&rcv_buffer, i_end - i_begin, MPI_DOUBLE, source_rank, sender_tag, MPI_COMM_WORLD, &current_request);
        int bi = 0;
        for (int i = i_begin; i < i_end; i++)
        {
            f(i,j_fixed) = rcv_buffer[bi];
            bi++;
        }
    }
		return current_request;
}
MPI_Request StaggeredGrid::receive_boundary_vertical_g(int sender_tag, int j_fixed, int i_begin, int i_end, int source_rank, bool do_nothing)
{

		MPI_Request current_request;
    if (!do_nothing)
    {
        std::vector<double> rcv_buffer;
        MPI_Irecv(&rcv_buffer, i_end - i_begin, MPI_DOUBLE, source_rank, sender_tag, MPI_COMM_WORLD, &current_request);
        int bi = 0;
        for (int i = i_begin; i < i_end; i++)
        {
            g(i,j_fixed) = rcv_buffer[bi];
            bi++;
        }
    }
		return current_request;
}
MPI_Request StaggeredGrid::receive_boundary_vertical_p(int sender_tag, int j_fixed, int i_begin, int i_end, int source_rank, bool do_nothing)
{

		MPI_Request current_request;
    if (!do_nothing)
    {
        std::vector<double> rcv_buffer;
        MPI_Irecv(&rcv_buffer, i_end - i_begin, MPI_DOUBLE, source_rank, sender_tag, MPI_COMM_WORLD, &current_request);
        int bi = 0;
        for (int i = i_begin; i < i_end; i++)
        {
            p(i,j_fixed) = rcv_buffer[bi];
            bi++;
        }
    }
		return current_request;
}

MPI_Request StaggeredGrid::receive_boundary_horizontal_u(int sender_tag, int i_fixed, int j_begin, int j_end, int source_rank, bool do_nothing)
{
		MPI_Request current_request;
    if (!do_nothing)
    {
        std::vector<double> rcv_buffer;
        MPI_Irecv(&rcv_buffer, j_end - j_begin, MPI_DOUBLE, source_rank, sender_tag, MPI_COMM_WORLD, &current_request);
        int bj = 0;
        for (int j = j_begin; j < j_end; j++)
        {
            u(i_fixed,j) = rcv_buffer[bj];
            bj++;
        }
    }
		return current_request;
}
MPI_Request StaggeredGrid::receive_boundary_horizontal_v(int sender_tag, int i_fixed, int j_begin, int j_end, int source_rank, bool do_nothing)
{
		MPI_Request current_request;
    if (!do_nothing)
    {
        std::vector<double> rcv_buffer;
        MPI_Irecv(&rcv_buffer, j_end - j_begin, MPI_DOUBLE, source_rank, sender_tag, MPI_COMM_WORLD, &current_request);
        int bj = 0;
        for (int j = j_begin; j < j_end; j++)
        {
            v(i_fixed,j) = rcv_buffer[bj];
            bj++;
        }
    }
		return current_request;
}
MPI_Request StaggeredGrid::receive_boundary_horizontal_f(int sender_tag, int i_fixed, int j_begin, int j_end, int source_rank, bool do_nothing)
{
		MPI_Request current_request;
    if (!do_nothing)
    {
        std::vector<double> rcv_buffer;
        MPI_Irecv(&rcv_buffer, j_end - j_begin, MPI_DOUBLE, source_rank, sender_tag, MPI_COMM_WORLD, &current_request);
        int bj = 0;
        for (int j = j_begin; j < j_end; j++)
        {
            f(i_fixed,j) = rcv_buffer[bj];
            bj++;
        }
    }
		return current_request;
}
MPI_Request StaggeredGrid::receive_boundary_horizontal_g(int sender_tag, int i_fixed, int j_begin, int j_end, int source_rank, bool do_nothing)
{
		MPI_Request current_request;
    if (!do_nothing)
    {
        std::vector<double> rcv_buffer;
        MPI_Irecv(&rcv_buffer, j_end - j_begin, MPI_DOUBLE, source_rank, sender_tag, MPI_COMM_WORLD, &current_request);
        int bj = 0;
        for (int j = j_begin; j < j_end; j++)
        {
            g(i_fixed,j) = rcv_buffer[bj];
            bj++;
        }
    }
		return current_request;
}
MPI_Request StaggeredGrid::receive_boundary_horizontal_p(int sender_tag, int i_fixed, int j_begin, int j_end, int source_rank, bool do_nothing)
{
		MPI_Request current_request;
    if (!do_nothing)
    {
        std::vector<double> rcv_buffer;
        MPI_Irecv(&rcv_buffer, j_end - j_begin, MPI_DOUBLE, source_rank, sender_tag, MPI_COMM_WORLD, &current_request);
        int bj = 0;
        for (int j = j_begin; j < j_end; j++)
        {
            p(i_fixed,j) = rcv_buffer[bj];
            bj++;
        }
    }
		return current_request;
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
