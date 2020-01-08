#include <precice/SolverInterface.hpp>
#include <assert.h>
#include <iostream>

#include <vector>

class Adapter
{
    private:
    precice::SolverInterface precice;
    int dim, meshID, vertexSize, readDataID, writeDataID, rank, size;
    int* vertexIDs;
    double dt, precice_dt;
    // double* coords;
    std::vector<double> coords;
    std::vector<int> vertexIDs2;
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
    Adapter(std::string participantName, std::string preciceConfigFile, int rank, int size, int vertexSize, std::string readDataName, std::string writeDataName)
    : participantName(participantName), rank(rank), size(size), vertexSize(vertexSize), readDataName(readDataName), writeDataName(writeDataName), precice(participantName, rank, size) // preciceConfigFile(preciceConfigFile)
    {
        precice.configure(preciceConfigFile);
    }

    void initialize(std::string participantMesh, std::vector<double> & coords)
    {
	    std::cout << "initialize: start" << std::endl;
        std::cout << precice.getDimensions() << std::endl;
        meshID = precice.getMeshID(participantMesh);
	    std::cout << meshID << std::endl;
        writeDataID = precice.getDataID(writeDataName, meshID);
	    std::cout << writeDataID << std::endl;
        readDataID = precice.getDataID(readDataName, meshID);
	    std::cout << readDataID << std::endl;
        std::cout << precice.getDimensions() << std::endl;
        dim = precice.getDimensions();
        // assert( dim == 2 );
        // coords = new double[vertexSize*2]; // input from settings
        // double* coords2 = &coords[];
        // vertexIDs = new int[vertexSize];
        int vertexIDs_[vertexSize] = {0};
        vertexIDs = vertexIDs;
        vertexIDs2 = std::vector<int>(vertexSize,0);
        std::cout << "vertexSize: " << vertexSize << " sizeof(vertexIDs): " << sizeof(vertexIDs)/sizeof(int) << "Coords: " << coords.size() << ", vs2: " << vertexIDs2.size() << std::endl;
        precice.setMeshVertices(meshID,vertexSize,coords.data(),vertexIDs2.data());
        std::cout << "vertexSize: " << vertexSize << " sizeof(vertexIDs): " << sizeof(vertexIDs)/sizeof(int) << std::endl;
        // delete[] coords; // wofür? by HENRIK
        // writeData = new double[vertexSize];
        // readData = new double[vertexSize];
        precice_dt = precice.initialize();

        std::cout << precice.getMeshVertexSize(meshID) << std::endl;
    }

    double get_dt(double dt_in)
    {
        dt = std::min(dt_in, precice_dt);
        return dt;
    }

    void readData(std::vector<double> & readData)
    {
		std::cout << "run: pre read, readData size:" << readData.size() << ", vertexSize: " << vertexSize << ", readDataID: " << readDataID << ", vertexIDs: " << sizeof(vertexIDs)/sizeof(int) <<  std::endl;
        precice.readBlockVectorData(readDataID, vertexSize, vertexIDs2.data(), readData.data());
    }

    void writeData(std::vector<double> & writeData)
    {
        precice.writeBlockVectorData(writeDataID, vertexSize, vertexIDs, writeData.data());
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
