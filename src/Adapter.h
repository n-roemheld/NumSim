#include <precice/SolverInterface.hpp>
#include <assert.h>
// #include <vector>

class Adapter
{
    private:
    precice::SolverInterface precice;
    int dim, meshID, vertexSize, readDataID, writeDataID, rank, size;
    int* vertexIDs;
    double dt, precice_dt;
    double* coords;
    // double* readData; // sync with staggered grid
    // double* writeData; // 
    std::string readDataName, writeDataName;
    std::string participantName, preciceConfigFile; // not needed?
    // ...

    // Shorthand definition ofr preCICE constants
	// Read cowid  (= co + w + i + d) as
	// constant "write initial data"
	const std::string& cowid = precice::constants::actionWriteInitialData();
	const std::string& coric = precice::constants::actionReadIterationCheckpoint();
	const std::string& cowic = precice::constants::actionWriteIterationCheckpoint();

    public:
    Adapter(std::string participantName, std::string preciceConfigFile, int rank, int size, int vertexSize, std::string readDataName, std::string writeDataName)
    : participantName(participantName), preciceConfigFile(preciceConfigFile), rank(rank), size(size), vertexSize(vertexSize), readDataName(readDataName), writeDataName(writeDataName), precice(participantName, rank, size)
    {
        precice.configure(preciceConfigFile);
    }

    void initialize(std::string participantMesh, double *coords)
    {
        dim = precice.getDimensions();
        assert( dim == 2 );
        meshID = precice.getMeshID(participantMesh);
        // coords = new double[vertexSize*dim]; // get from settings
        vertexIDs = new int[vertexSize];
        precice.setMeshVertices(meshID,vertexSize,coords,vertexIDs);
        delete[] coords;
        writeDataID = precice.getDataID(writeDataName, meshID); 
        readDataID = precice.getDataID(readDataName, meshID); 
        // writeData = new double[vertexSize];
        // readData = new double[vertexSize];
        precice_dt = precice.initialize();
    }

    double get_dt(double dt_in)
    {
        dt = std::min(dt_in, precice_dt);
        return dt;
    }

    void readData(double *readData)
    {
        precice.readBlockVectorData(readDataID, vertexSize, vertexIDs, readData);
    }

    void writeData(double *writeData)
    {
        precice.writeBlockVectorData(writeDataID, vertexSize, vertexIDs, writeData);
    }

    void advance()
    {
        precice_dt = precice.advance(dt);
    }

    void finalize()
    {
        precice.finalize();
    }

    int getVertexSize()
    {
        return vertexSize;
    }
};