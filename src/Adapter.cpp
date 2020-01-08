#include "Adapter.h"

Adapter::Adapter(std::string participantName, std::string preciceConfigFile, int rank, int size, int vertexSize, std::string readDataName, std::string writeDataName)
: participantName(participantName), rank(rank), size(size), vertexSize(vertexSize), readDataName(readDataName), writeDataName(writeDataName), precice(participantName, rank, size) // preciceConfigFile(preciceConfigFile)
{
  std::cout << "Adapter Constructor" << std::endl;
    precice.configure(preciceConfigFile);
    std::cout << precice.getDimensions() << std::endl;
}

void Adapter::initialize(std::string participantMesh, std::vector<double> & coords)
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
    vertexIDs = new int[vertexSize];
    precice.setMeshVertices(meshID,vertexSize,coords.data(),vertexIDs);
    // delete[] coords; // wofür? by HENRIK
    // writeData = new double[vertexSize];
    // readData = new double[vertexSize];
    precice_dt = precice.initialize();
}

double Adapter::get_dt(double dt_in)
{
    dt = std::min(dt_in, precice_dt);
    return dt;
}

void Adapter::readData(std::vector<double> & readData)
{
    precice.readBlockVectorData(readDataID, vertexSize, vertexIDs, readData.data());
}

void Adapter::writeData(std::vector<double> & writeData)
{
    precice.writeBlockVectorData(writeDataID, vertexSize, vertexIDs, writeData.data());
}

void Adapter::advance()
{
    precice_dt = precice.advance(dt);
}

void Adapter::finalize()
{
    precice.finalize();
}

int Adapter::getVertexSize()
{
    return vertexSize;
}
