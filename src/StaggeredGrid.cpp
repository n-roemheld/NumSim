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
    return partitioning_.is_boundary(direction) == 1;
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


// void StaggeredGrid::receive_boundary(int sender_tag, int index_fixed, int index_begin, int index_end, int source_rank, bool do_nothing, char function, bool horizontal_communication)
// {
//     if (!do_nothing)
//     {
//         std::vector<double> rcv_buffer(index_end - index_begin);
//         MPI_Recv(&rcv_buffer, index_end - index_begin, MPI_DOUBLE, source_rank, sender_tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//         write_recv_buffer(index_fixed, index_begin, index_end, function, horizontal_communication, rcv_buffer);
//     }
// };


void StaggeredGrid::write_recv_buffer(int index_fixed, int index_begin, int index_end, char function, bool horizontal_communication, std::vector<double> rcv_buffer)
{
    int buffer_index = 0;
    int index1;
    int index2;
    for (int index = index_begin; index < index_end; index++)
    {
        if (horizontal_communication) {
            index1 = index_fixed;
            index2 = index;
        }
        if (function == 'p') {
            p(index1,index2) = rcv_buffer.at(buffer_index);
        } else if (function == 'u') {
            u(index1,index2) = rcv_buffer.at(buffer_index);
        } else if (function == 'v') {
            v(index1,index2) = rcv_buffer.at(buffer_index);
        } else if (function == 'f') {
            f(index1,index2) = rcv_buffer.at(buffer_index);
        } else if (function == 'g') {
            g(index1,index2) = rcv_buffer.at(buffer_index);
        }
        buffer_index++;
    }
    
};

// void StaggeredGrid::send_boundary(int receiver_tag, int index_fixed, int index_begin, int index_end, int target_rank, bool do_nothing, char function, bool horizontal_communication)
// {
//     if (!do_nothing)
//     {
//         std::vector<double> send_buffer(index_end - index_begin);
//         send_buffer = get_send_buffer(index_fixed, index_begin, index_end, function, horizontal_communication);
//         MPI_Send(&send_buffer, index_end - index_begin, MPI_DOUBLE, target_rank, receiver_tag, MPI_COMM_WORLD);
//     }
// };

std::vector<double> StaggeredGrid::get_send_buffer(int index_fixed, int index_begin, int index_end, char function, bool horizontal_communication)
{
    std::vector<double> send_buffer(index_end - index_begin);
    int buffer_index = 0;
    int index1;
    int index2;
    for (int index = index_begin; index < index_end; index++)
    {
        if (horizontal_communication) {
            index1 = index_fixed;
            index2 = index;
        }
        if (function == 'p') {
            send_buffer.at(buffer_index) = p(index1,index2);
        } else if (function == 'u') {
            send_buffer.at(buffer_index) = u(index1,index2);
        } else if (function == 'v') {
            send_buffer.at(buffer_index) = v(index1,index2);
        } else if (function == 'f') {
            send_buffer.at(buffer_index) = f(index1,index2);
        } else if (function == 'g') {
            send_buffer.at(buffer_index) = g(index1,index2);
        } else {
            std::cout << "Function not recognized!" << std::endl;
        }
        buffer_index++;
    }
    return send_buffer;
}

void StaggeredGrid::velocity_horizontal_communication(char u_or_f, char v_or_g)
{
    
    // neigbor indices
    int right = 1;
    int left = 3;

    // MPI tags of information sent (communication directions)
    int to_right = 1;
    int to_left = 3;

    // Index bounds of the values to send
    int uJ_start = uJBegin();
    int uJ_end = uJEnd();

    int vJ_start = vJBegin();
    int vJ_end = uJEnd()-1;

    std::vector<MPI_Request> requests;

    // Communicate to right side (u and v, send and receive)
    if (!is_boundary(right)){
        // send u
        std::vector<double> send_buffer_u = get_send_buffer(uIEnd()-2, uJ_start, uJ_end, u_or_f, true);
        int nValues = send_buffer_u.size();
        requests.emplace_back();
        MPI_Isend(send_buffer_u.data(), nValues, MPI_DOUBLE, rank_neighbor(right), to_right, MPI_COMM_WORLD, &requests.back());
        // receive u
        std::vector<double> recv_buffer_u(nValues);
        requests.emplace_back();
        MPI_Irecv(recv_buffer_u.data(), nValues, MPI_DOUBLE, rank_neighbor(right), to_left, MPI_COMM_WORLD, &requests.back());
        write_recv_buffer(uIEnd()-2, uJ_start, uJ_end, u_or_f, true, recv_buffer_u);

        // send v
        std::vector<double> send_buffer_v = get_send_buffer(vIEnd()-1, vJ_start, vJ_end, v_or_g, true);
        nValues = send_buffer_v.size();
        requests.emplace_back();
        MPI_Isend(send_buffer_v.data(), nValues, MPI_DOUBLE, rank_neighbor(right), 4 + to_right, MPI_COMM_WORLD, &requests.back());
        // receive v
        std::vector<double> recv_buffer_v(nValues);
        requests.emplace_back();
        MPI_Irecv(recv_buffer_v.data(), nValues, MPI_DOUBLE, rank_neighbor(right), 4 + to_left, MPI_COMM_WORLD, &requests.back());
        write_recv_buffer(vIEnd(), vJ_start, vJ_end, v_or_g, true, recv_buffer_v);
    }   
    // Communicate to left side (u and v, send and receive)
    if (!is_boundary(left)) {
        // send u
        std::vector<double> send_buffer_u = get_send_buffer(uIBegin(), uJ_start, uJ_end, u_or_f, true);
        int nValues = send_buffer_u.size();
        requests.emplace_back();
        MPI_Isend(send_buffer_u.data(), nValues, MPI_DOUBLE, rank_neighbor(left), to_left, MPI_COMM_WORLD, &requests.back());
        // receive u
        std::vector<double> recv_buffer_u(nValues);
        requests.emplace_back();
        MPI_Irecv(recv_buffer_u.data(), nValues, MPI_DOUBLE, rank_neighbor(left), to_right, MPI_COMM_WORLD, &requests.back());
        write_recv_buffer(uIBegin()-1, uJ_start, uJ_end, u_or_f, true, recv_buffer_u);

        // send v
        std::vector<double> send_buffer_v = get_send_buffer(vIBegin(), vJ_start, vJ_end, v_or_g, true);
        nValues = send_buffer_v.size();
        requests.emplace_back();
        MPI_Isend(send_buffer_v.data(), nValues, MPI_DOUBLE, rank_neighbor(left), 4 + to_left, MPI_COMM_WORLD, &requests.back());
        // receive v
        std::vector<double> recv_buffer_v(nValues);
        requests.emplace_back();
        MPI_Irecv(recv_buffer_v.data(), nValues, MPI_DOUBLE, rank_neighbor(left), 4 + to_right, MPI_COMM_WORLD, &requests.back());
        write_recv_buffer(vIBegin()-1, vJ_start, vJ_end, v_or_g, true, recv_buffer_v);
    }

    MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);
}

void StaggeredGrid::velocity_vertical_communication(char u_or_f, char v_or_g)
{
    
    // neigbor indices
    int below = 0;
    int above = 2;

    // MPI tags of information sent (communication directions)
    int down = 0;
    int up = 2;

    // Index bounds of the values to send
    int uI_start = uIBegin();
    int uI_end = uIEnd();

    int vI_start = vIBegin()-1;
    int vI_end = uIEnd();

    std::vector<MPI_Request> requests;

    // Communicate to below side (u and v, send and receive)
    if (!is_boundary(below)) {
        // send u
        std::vector<double> send_buffer_u = get_send_buffer(uJBegin(), uI_start, uI_end, u_or_f, false);
        int nValues = send_buffer_u.size();
        requests.emplace_back();
        MPI_Isend(send_buffer_u.data(), nValues, MPI_DOUBLE, rank_neighbor(below), down, MPI_COMM_WORLD, &requests.back());
        // receive u
        std::vector<double> recv_buffer_u(nValues);
        requests.emplace_back();
        MPI_Irecv(recv_buffer_u.data(), nValues, MPI_DOUBLE, rank_neighbor(below), up, MPI_COMM_WORLD, &requests.back());
        write_recv_buffer(uJBegin()-1, uI_start, uI_end, u_or_f, false, recv_buffer_u);

        // send v
        std::vector<double> send_buffer_v = get_send_buffer(vJBegin(), vI_start, vI_end, v_or_g, false);
        nValues = send_buffer_v.size();
        requests.emplace_back();
        MPI_Isend(send_buffer_v.data(), nValues, MPI_DOUBLE, rank_neighbor(below), 4 + down, MPI_COMM_WORLD, &requests.back());
        // receive v
        std::vector<double> recv_buffer_v(nValues);
        requests.emplace_back();
        MPI_Irecv(recv_buffer_v.data(), nValues, MPI_DOUBLE, rank_neighbor(below), 4 + up, MPI_COMM_WORLD, &requests.back());
        write_recv_buffer(vJBegin()-1, vI_start, vI_end, v_or_g, false, recv_buffer_v);
    }   
    // Communicate to above side (u and v, send and receive)
    if (!is_boundary(above)) {
        // send u
        std::vector<double> send_buffer_u = get_send_buffer(uJEnd(), uI_start, uI_end, u_or_f, false);
        int nValues = send_buffer_u.size();
        requests.emplace_back();
        MPI_Isend(send_buffer_u.data(), nValues, MPI_DOUBLE, rank_neighbor(above), up, MPI_COMM_WORLD, &requests.back());
        // receive u
        std::vector<double> recv_buffer_u(nValues);
        requests.emplace_back();
        MPI_Irecv(recv_buffer_u.data(), nValues, MPI_DOUBLE, rank_neighbor(above), down, MPI_COMM_WORLD, &requests.back());
        write_recv_buffer(uJEnd(), uI_start, uI_end, u_or_f, false, recv_buffer_u);

        // send v
        std::vector<double> send_buffer_v = get_send_buffer(vJEnd()-2, vI_start, vI_end, v_or_g, false);
        nValues = send_buffer_v.size();
        requests.emplace_back();
        MPI_Isend(send_buffer_v.data(), nValues, MPI_DOUBLE, rank_neighbor(above), 4 + up, MPI_COMM_WORLD, &requests.back());
        // receive v
        std::vector<double> recv_buffer_v(nValues);
        requests.emplace_back();
        MPI_Irecv(recv_buffer_v.data(), nValues, MPI_DOUBLE, rank_neighbor(above), 4 + down, MPI_COMM_WORLD, &requests.back());
        write_recv_buffer(vJEnd()-1, vI_start, vI_end, v_or_g, false, recv_buffer_v);
    }

    MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);
}



void StaggeredGrid::pressure_communication()
{
    // neigbor indices
    int below = 0;
    int above = 2;
    int right = 1;
    int left = 3;

    // MPI tags of information sent (communication directions)
    int down = 0;
    int to_right = 1;
    int up = 2;
    int to_left = 3;

    std::vector<MPI_Request> requests;

    if (!is_boundary(right)) {
        // send
        std::vector<double> send_buffer = get_send_buffer(pIEnd()-1, pJBegin(), pJEnd(), 'p', false);
        int nValues = send_buffer.size();
        requests.emplace_back();
        MPI_Isend(send_buffer.data(), nValues, MPI_DOUBLE, rank_neighbor(right), to_right, MPI_COMM_WORLD, &requests.back());
        // receive 
        std::vector<double> recv_buffer(nValues);
        requests.emplace_back();
        MPI_Irecv(recv_buffer.data(), nValues, MPI_DOUBLE, rank_neighbor(right), to_left, MPI_COMM_WORLD, &requests.back());
        write_recv_buffer(pIEnd(), pJBegin(), pJEnd(), 'p', false, recv_buffer);   
        
    }
    if (!is_boundary(left)) {
        // send
        std::vector<double> send_buffer = get_send_buffer(pIBegin(), pJBegin(), pJEnd(), 'p', false);
        int nValues = send_buffer.size();
        requests.emplace_back();
        MPI_Isend(send_buffer.data(), nValues, MPI_DOUBLE, rank_neighbor(left), to_left, MPI_COMM_WORLD, &requests.back());
        // receive 
        std::vector<double> recv_buffer(nValues);
        requests.emplace_back();
        MPI_Irecv(recv_buffer.data(), nValues, MPI_DOUBLE, rank_neighbor(left), to_right, MPI_COMM_WORLD, &requests.back());
        write_recv_buffer(pIBegin()-1, pJBegin(), pJEnd(), 'p', false, recv_buffer);   
        
    }
    if (!is_boundary(below)) {
        // send
        std::vector<double> send_buffer = get_send_buffer(pJBegin(), pIBegin(), pIEnd(), 'p', false);
        int nValues = send_buffer.size();
        requests.emplace_back();
        MPI_Isend(send_buffer.data(), nValues, MPI_DOUBLE, rank_neighbor(below), down, MPI_COMM_WORLD, &requests.back());
        // receive 
        std::vector<double> recv_buffer(nValues);
        requests.emplace_back();
        MPI_Irecv(recv_buffer.data(), nValues, MPI_DOUBLE, rank_neighbor(below), up, MPI_COMM_WORLD, &requests.back());
        write_recv_buffer(pJBegin()-1, pIBegin(), pIEnd(), 'p', false, recv_buffer);   
        
    }
    if (!is_boundary(above)) {
        // send
        std::vector<double> send_buffer = get_send_buffer(pJEnd()-1, pIBegin(), pIEnd(), 'p', false);
        int nValues = send_buffer.size();
        requests.emplace_back();
        MPI_Isend(send_buffer.data(), nValues, MPI_DOUBLE, rank_neighbor(above), up, MPI_COMM_WORLD, &requests.back());
        // receive 
        std::vector<double> recv_buffer(nValues);
        requests.emplace_back();
        MPI_Irecv(recv_buffer.data(), nValues, MPI_DOUBLE, rank_neighbor(above), down, MPI_COMM_WORLD, &requests.back());
        write_recv_buffer(pJEnd(), pIBegin(), pIEnd(), 'p', false, recv_buffer);       
    }    
    MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);
};