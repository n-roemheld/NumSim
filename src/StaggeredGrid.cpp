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


// vector<double> StaggeredGrid::exchangeValues(vector<double> values, int receiverRank, vector<MPI_Request> &requests) {
//     int nValues = values.size();
//     std::vector<double> receiveBuffer(nValues);
//     requests.emplace_back();
//     MPI_Isend(values.data(), nValues, MPI_DOUBLE, receiverRank, 0,
//               MPI_COMM_WORLD, &requests.back());
//     requests.emplace_back();
//     MPI_Irecv(receiveBuffer.data(), nValues, MPI_DOUBLE, receiverRank, 0,
//               MPI_COMM_WORLD, &requests.back());
//     return receiveBuffer;
// }
//
// void Communication::communicate(FieldVariable variable, std::string type) {
//     std::vector<MPI_Request> requests;
//     std::array<std::vector<double>, 4> sendBuffers;
//     std::array<std::vector<double>, 4> receiveBuffers;
//     // 0 is top
//     // 1 is bottom
//     // 2 is left
//     // 3 is right
//     sendBuffers[0] = std::vector<double>(variable.size()[0]);
//     sendBuffers[1] = std::vector<double>(variable.size()[0]);
//     sendBuffers[2] = std::vector<double>(variable.size()[1]);
//     sendBuffers[3] = std::vector<double>(variable.size()[1]);
//     if (type == "p") {
//         for (int i = 0; i < sendBuffers[0].size(); i++) {
//             sendBuffers[0][i] = variable.operator()(i,
//                                                     discretization_.get()->pJEnd());
//             sendBuffers[1][i] = variable.operator()(i,
//                                                     discretization_.get()->pJBegin());
//         }
//         for (int j = 0; j < sendBuffers[2].size(); j++) {
//             sendBuffers[2][j] = variable.operator()(discretization_.get()->pIBegin(),
//                                                     j);
//             sendBuffers[3][j] = variable.operator()(discretization_.get()->pIEnd(),
//                                                     j);
//         }
//     } else if (type == "u" or type == "f") {
//         for (int i = 0; i < sendBuffers[0].size(); i++) {
//             sendBuffers[0][i] = variable.operator()(i,
//                                                     discretization_.get()->uJEnd());
//             sendBuffers[1][i] = variable.operator()(i,
//                                                     discretization_.get()->uJBegin());
//         }
//         for (int j = 0; j < sendBuffers[2].size(); j++) {
//             sendBuffers[2][j] = variable.operator()(discretization_.get()->uIBegin(),
//                                                     j);
//             sendBuffers[3][j] = variable.operator()(discretization_.get()->uIEnd(),
//                                                     j);
//         }
//     } else if (type == "v" or type == "g") {
//         for (int i = 0; i < sendBuffers[0].size(); i++) {
//             sendBuffers[0][i] = variable.operator()(i,
//                                                     discretization_.get()->vJEnd());
//             sendBuffers[1][i] = variable.operator()(i,
//                                                     discretization_.get()->vJBegin());
//         }
//         for (int j = 0; j < sendBuffers[2].size(); j++) {
//             sendBuffers[2][j] = variable.operator()(discretization_.get()->vIBegin(),
//                                                     j);
//             sendBuffers[3][j] = variable.operator()(discretization_.get()->vIEnd(),
//                                                     j);
//         }
//     }
//     if (partitioning_.get()->getRankOfTopNeighbour() != -1) {
//         receiveBuffers[0] = exchangeValues(sendBuffers[0], partitioning_.get()->getRankOfTopNeighbour(), requests);
//     }
//     if (partitioning_.get()->getRankOfBottomNeighbour() != -1) {
//         receiveBuffers[1] = exchangeValues(sendBuffers[1], partitioning_.get()->getRankOfBottomNeighbour(), requests);
//     }
//     if (partitioning_.get()->getRankOfLeftNeighbour() != -1) {
//         receiveBuffers[2] = exchangeValues(sendBuffers[2], partitioning_.get()->getRankOfLeftNeighbour(), requests);
//     }
//     if (partitioning_.get()->getRankOfRightNeighbour() != -1) {
//         receiveBuffers[3] = exchangeValues(sendBuffers[3], partitioning_.get()->getRankOfRightNeighbour(), requests);
//     }
//     MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);
//
//     // Writing back top border
//     if (type == "p") {
//         if (partitioning_.get()->getRankOfTopNeighbour() != -1) {
//             for (int i = 0; i < receiveBuffers[0].size(); i++) {
//                 discretization_.get()->p(i,
//                                          discretization_.get()->pJEnd() + 1) = receiveBuffers[0][i];
//             }
//         }
//         if (partitioning_.get()->getRankOfBottomNeighbour() != -1) {
//             for (int i = 0; i < receiveBuffers[1].size(); i++) {
//                 discretization_.get()->p(i,
//                                          discretization_.get()->pJBegin() - 1) = receiveBuffers[1][i];
//             }
//         }
//         if (partitioning_.get()->getRankOfLeftNeighbour() != -1) {
//             for (int j = 0; j < receiveBuffers[2].size(); j++) {
//                 discretization_.get()->p(discretization_.get()->pIBegin() - 1, j) = receiveBuffers[2][j];
//             }
//         }
//         if (partitioning_.get()->getRankOfRightNeighbour() != -1) {
//             for (int j = 0; j < receiveBuffers[3].size(); j++) {
//                 discretization_.get()->p(discretization_.get()->pIEnd() + 1, j) = receiveBuffers[3][j];
//             }
//         }
//
//     } else if (type == "u") {
//         if (partitioning_.get()->getRankOfTopNeighbour() != -1) {
//             for (int i = 0; i < receiveBuffers[0].size(); i++) {
//                 discretization_.get()->u(i,
//                                          discretization_.get()->uJEnd() + 1) = receiveBuffers[0][i];
//             }
//         }
//         if (partitioning_.get()->getRankOfBottomNeighbour() != -1) {
//             for (int i = 0; i < receiveBuffers[1].size(); i++) {
//                 discretization_.get()->u(i,
//                                          discretization_.get()->uJBegin() - 1) = receiveBuffers[1][i];
//             }
//         }
//         if (partitioning_.get()->getRankOfLeftNeighbour() != -1) {
//             for (int j = 0; j < receiveBuffers[2].size(); j++) {
//                 discretization_.get()->u(discretization_.get()->uIBegin() - 1, j) = receiveBuffers[2][j];
//             }
//         }
//         if (partitioning_.get()->getRankOfRightNeighbour() != -1) {
//             for (int j = 0; j < receiveBuffers[3].size(); j++) {
//                 discretization_.get()->u(discretization_.get()->uIEnd() + 1, j) = receiveBuffers[3][j];
//             }
//         }
//     } else if (type == "f") {
//         if (partitioning_.get()->getRankOfTopNeighbour() != -1) {
//             for (int i = 0; i < receiveBuffers[0].size(); i++) {
//                 discretization_.get()->f(i,
//                                          discretization_.get()->uJEnd() + 1) = receiveBuffers[0][i];
//             }
//         }
//         if (partitioning_.get()->getRankOfBottomNeighbour() != -1) {
//             for (int i = 0; i < receiveBuffers[1].size(); i++) {
//                 discretization_.get()->f(i,
//                                          discretization_.get()->uJBegin() - 1) = receiveBuffers[1][i];
//             }
//         }
//         if (partitioning_.get()->getRankOfLeftNeighbour() != -1) {
//             for (int j = 0; j < receiveBuffers[2].size(); j++) {
//                 discretization_.get()->f(discretization_.get()->uIBegin() - 1, j) = receiveBuffers[2][j];
//             }
//         }
//         if (partitioning_.get()->getRankOfRightNeighbour() != -1) {
//             for (int j = 0; j < receiveBuffers[3].size(); j++) {
//                 discretization_.get()->f(discretization_.get()->uIEnd() + 1, j) = receiveBuffers[3][j];
//             }
//         }
//     } else if (type == "v") {
//         if (partitioning_.get()->getRankOfTopNeighbour() != -1) {
//             for (int i = 0; i < receiveBuffers[0].size(); i++) {
//                 discretization_.get()->v(i,
//                                          discretization_.get()->vJEnd() + 1) = receiveBuffers[0][i];
//             }
//         }
//         if (partitioning_.get()->getRankOfBottomNeighbour() != -1) {
//             for (int i = 0; i < receiveBuffers[1].size(); i++) {
//                 discretization_.get()->v(i,
//                                          discretization_.get()->vJBegin() - 1) = receiveBuffers[1][i];
//             }
//         }
//         if (partitioning_.get()->getRankOfLeftNeighbour() != -1) {
//             for (int j = 0; j < receiveBuffers[2].size(); j++) {
//                 discretization_.get()->v(discretization_.get()->vIBegin() - 1, j) = receiveBuffers[2][j];
//             }
//         }
//         if (partitioning_.get()->getRankOfRightNeighbour() != -1) {
//             for (int j = 0; j < receiveBuffers[3].size(); j++) {
//                 discretization_.get()->v(discretization_.get()->vIEnd() + 1,
//                                          j) = receiveBuffers[3][j];
//             }
//         }
//
//     } else if (type == "g") {
//         if (partitioning_.get()->getRankOfTopNeighbour() != -1) {
//             for (int i = 0; i < receiveBuffers[0].size(); i++) {
//                 discretization_.get()->g(i,
//                                          discretization_.get()->vJEnd() + 1) = receiveBuffers[0][i];
//             }
//         }
//         if (partitioning_.get()->getRankOfBottomNeighbour() != -1) {
//             for (int i = 0; i < receiveBuffers[1].size(); i++) {
//                 discretization_.get()->g(i,
//                                          discretization_.get()->vJBegin() - 1) = receiveBuffers[1][i];
//             }
//         }
//         if (partitioning_.get()->getRankOfLeftNeighbour() != -1) {
//             for (int j = 0; j < receiveBuffers[2].size(); j++) {
//                 discretization_.get()->g(discretization_.get()->vIBegin() - 1, j) = receiveBuffers[2][j];
//             }
//         }
//         if (partitioning_.get()->getRankOfRightNeighbour() != -1) {
//             for (int j = 0; j < receiveBuffers[3].size(); j++) {
//                 discretization_.get()->g(discretization_.get()->vIEnd() + 1,
//                                          j) = receiveBuffers[3][j];
//             }
//         }
//
//     }
// 	}

	////////////////////////////////////////////////////////////////////////

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
				else
				{
					index1 = index;
					index2 = index_fixed;
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
				else
				{
					index1 = index;
					index2 = index_fixed;
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

		// unnecessary stuff

		int buffer= 0;
		requests.emplace_back();
		MPI_Isend(&buffer, 1, MPI_INT, 0, 7, MPI_COMM_WORLD, &requests.back());

		requests.emplace_back();
		int buffer2;
		MPI_Irecv(&buffer2, 1, MPI_INT, 0, 7, MPI_COMM_WORLD, &requests.back());


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
                                        write_recv_buffer(uIEnd()-1, uJ_start, uJ_end, u_or_f, true, recv_buffer_u);

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
		int MPI_rank;
		MPI_Comm_rank(MPI_COMM_WORLD, &MPI_rank);

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

		// unnecessary stuff

	  int buffer= 0;
	  requests.emplace_back();
	  MPI_Isend(&buffer, 1, MPI_INT, 0, 7, MPI_COMM_WORLD, &requests.back());

	  requests.emplace_back();
	  int buffer2;
	  MPI_Irecv(&buffer2, 1, MPI_INT, 0, 7, MPI_COMM_WORLD, &requests.back());


    // Communicate to below side (u and v, send and receive)
    if (!is_boundary(below)) {
        // send u
        std::vector<double> send_buffer_u = get_send_buffer(uJBegin(), uI_start, uI_end, u_or_f, false);
        int nValues = send_buffer_u.size();
        requests.emplace_back();
        MPI_Isend(send_buffer_u.data(), nValues, MPI_DOUBLE, rank_neighbor(below), down, MPI_COMM_WORLD, &requests.back());
				std::cout << MPI_rank << "send u"<< "tag: " << down << "size :" << nValues << std::endl;
        // receive u
        std::vector<double> recv_buffer_u(nValues);
        requests.emplace_back();
				std::cout<< MPI_rank << "receive u"<< "tag: " << up << "size :" << nValues << std::endl;
        MPI_Irecv(recv_buffer_u.data(), nValues, MPI_DOUBLE, rank_neighbor(below), up, MPI_COMM_WORLD, &requests.back());
                                        write_recv_buffer(uJBegin()-1, uI_start, uI_end, u_or_f, false, recv_buffer_u);

        // send v
        std::vector<double> send_buffer_v = get_send_buffer(vJBegin(), vI_start, vI_end, v_or_g, false);
        nValues = send_buffer_v.size();
        requests.emplace_back();
        MPI_Isend(send_buffer_v.data(), nValues, MPI_DOUBLE, rank_neighbor(below), 4 + down, MPI_COMM_WORLD, &requests.back());
				std::cout << MPI_rank << "send v"<< "tag: " << 4+down << "size :" << nValues << std::endl;
        // receive v
        std::vector<double> recv_buffer_v(nValues);
        requests.emplace_back();
				std::cout << MPI_rank << "receive v"<< "tag: " << 4+up << "size :" << nValues << std::endl;
        MPI_Irecv(recv_buffer_v.data(), nValues, MPI_DOUBLE, rank_neighbor(below), 4 + up, MPI_COMM_WORLD, &requests.back());
                                        write_recv_buffer(vJBegin()-1, vI_start, vI_end, v_or_g, false, recv_buffer_v);
    }
    // Communicate to above side (u and v, send and receive)
    if (!is_boundary(above)) {
        // send u
        std::vector<double> send_buffer_u = get_send_buffer(uJEnd()-1, uI_start, uI_end, u_or_f, false);
        int nValues = send_buffer_u.size();
        requests.emplace_back();
				std::cout << MPI_rank << "send u"<< "tag: " << up << "size :" << nValues << std::endl;
        MPI_Isend(send_buffer_u.data(), nValues, MPI_DOUBLE, rank_neighbor(above), up, MPI_COMM_WORLD, &requests.back());
        // receive u
        std::vector<double> recv_buffer_u(nValues);
        requests.emplace_back();
				std::cout << MPI_rank << "receive u"<< "tag: " << down << "size :" << nValues << std::endl;
        MPI_Irecv(recv_buffer_u.data(), nValues, MPI_DOUBLE, rank_neighbor(above), down, MPI_COMM_WORLD, &requests.back());
                                            write_recv_buffer(uJEnd(), uI_start, uI_end, u_or_f, false, recv_buffer_u);

        // send v
        std::vector<double> send_buffer_v = get_send_buffer(vJEnd()-2, vI_start, vI_end, v_or_g, false);
        nValues = send_buffer_v.size();
        requests.emplace_back();
				std::cout << MPI_rank << "send v"<< "tag: " << 4+up << "size :" << nValues << std::endl;
        MPI_Isend(send_buffer_v.data(), nValues, MPI_DOUBLE, rank_neighbor(above), 4 + up, MPI_COMM_WORLD, &requests.back());
        // receive v
        std::vector<double> recv_buffer_v(nValues);
        requests.emplace_back();
				std::cout << MPI_rank << "receive v"<< "tag: " << 4+down << "size :" << nValues << std::endl;
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


		// unnecessary stuff

	  int buffer= 0;
	  requests.emplace_back();
	  MPI_Isend(&buffer, 1, MPI_INT, 0, 7, MPI_COMM_WORLD, &requests.back());

	  requests.emplace_back();
	  int buffer2;
	  MPI_Irecv(&buffer2, 1, MPI_INT, 0, 7, MPI_COMM_WORLD, &requests.back());


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
