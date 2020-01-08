#include <precice/SolverInterface.hpp>
#include <assert.h>
#include <iostream>

#include <vector>

class Adapter
{
    private:
    precice::SolverInterface precice;
    int dim, meshID, vertexSize, readDataID, writeDataID, rank, size;
    std::vector<int> vertexIDs;
    double dt, precice_dt;
    // double* coords;
    std::vector<double> coords;
    // double* readData; // sync with staggered grid
    // double* writeData; //
    std::string readDataName, writeDataName;
    std::string participantName;
    // std::string preciceConfigFile; // not needed?
    // ...

    // Shorthand definition for preCICE constants
	// Read cowid  (= co + w + i + d) as
	// constant "write initial data"
	const std::string& cowid = precice::constants::actionWriteInitialData();
	const std::string& coric = precice::constants::actionReadIterationCheckpoint();
	const std::string& cowic = precice::constants::actionWriteIterationCheckpoint();

    public:
    Adapter(std::string participantName, std::string preciceConfigFile, int rank, int size, int vertexSize, std::string readDataName, std::string writeDataName);

    void initialize(std::string participantMesh, std::vector<double> coords);
    double get_dt(double dt_in);
    void readData(std::vector<double> readData);
    // void readData();

    void writeData(std::vector<double> writeData);
    void advance();
    void finalize();
    int getVertexSize();
};
