#include "settings.h"
#include <fstream>   // for file operations
#include <iostream>  // for cout

void Settings::loadFromFile(std::string filename) {
	// input file
	// open file
	std::ifstream file(filename.c_str(), std::ios::in);

	// check if file is open
	if (!file.is_open()) {
		std::cout << "Could not open parameter file \"" << filename << "\"."
				<< std::endl;
		return;
	}

	// loop over lines of file
	for (int lineNo = 0;; lineNo++) {
		// read line
		std::string line;
		getline(file, line);
//    while(line[0] == ' '  || line[0]== '\t')
//    {
//    	line.erase(0,1);
//    }
		// erase blank spaces
		if (line.find_first_of(" \t") != std::string::npos) {
			line.erase(0, line.find_first_not_of(" \t"));
		}

		if (line[0] != '#' && line.find_first_of('=') != std::string::npos)
		{
			std::string parameterName = line.substr(0, line.find_first_of('=')); //-1 ??
			// delete spaces after parameterName
			if (parameterName.find_first_of(" \t") != std::string::npos)
			{
				parameterName.erase(parameterName.find_first_of(" \t"));
			}

			std::string parameterValue = line.substr(
					line.find_first_of('=') + 1);
			// delete spaces before parameterValue
			if (parameterValue.find_first_of(" \t") != std::string::npos)
			{
				parameterValue.erase(0,
						parameterValue.find_first_not_of(" \t"));
			}
			// delete commands
			if (parameterValue.find_first_of("#") != std::string::npos)
			{
				parameterValue.erase(parameterValue.find_first_of("#"));
			}

			// delete spaces after parameterValue
			if (parameterValue.find_first_of(" \t") != std::string::npos)
			{
				parameterValue.erase(parameterValue.find_first_of(" \t"));
			}

			if (parameterName == "physicalSizeX")
			{
				physicalSize[0] = atof(parameterValue.c_str());
			}
			else if (parameterName == "physicalSizeY")
			{
				physicalSize[1] = atof(parameterValue.c_str());
			}
			else if (parameterName == "endTime")
			{
				endTime = atof(parameterValue.c_str());
			}
			else if (parameterName == "re") {
				re = int(atof(parameterValue.c_str()));
			}
			else if (parameterName == "prandtl")
			{
				prandtl = atof(parameterValue.c_str());
			}
			else if (parameterName == "beta")
			{
				beta = atof(parameterValue.c_str());
			}
			else if (parameterName == "gX")
			{
				g[0] = atof(parameterValue.c_str());
			}
			else if (parameterName == "gY")
			{
				g[1] = atof(parameterValue.c_str());
			}
			else if (parameterName == "dirichletBottomX")
			{
				dirichletBcBottom[0] = atof(parameterValue.c_str());
			}
			else if (parameterName == "dirichletBottomY")
			{
				dirichletBcBottom[1] = atof(parameterValue.c_str());
			}
			else if (parameterName == "dirichletTopX")
			{
				dirichletBcTop[0] = atof(parameterValue.c_str());
			}
			else if (parameterName == "dirichletTopY")
			{
				dirichletBcTop[1] = atof(parameterValue.c_str());
			}
			else if (parameterName == "dirichletLeftX")
			{
				dirichletBcLeft[0] = atof(parameterValue.c_str());
			}
			else if (parameterName == "dirichletLeftY")
			{
				dirichletBcLeft[1] = atof(parameterValue.c_str());
			}
			else if (parameterName == "dirichletRightX")
			{
				dirichletBcRight[0] = atof(parameterValue.c_str());
			}
			else if (parameterName == "dirichletRightY")
			{
				dirichletBcRight[1] = atof(parameterValue.c_str());
			}
			else if (parameterName == "nCellsX")
			{
				nCells[0] = int(atof(parameterValue.c_str()));
			}
			else if (parameterName == "nCellsY")
			{
				nCells[1] = int(atof(parameterValue.c_str()));
			}
			else if (parameterName == "useDonorCell")
			{
				if (parameterValue == "true")
				{
					useDonorCell = true;
				}
				else if (parameterValue == "false")
				{
					useDonorCell = false;
				}
				else
				{
					std::cout << "unknown useDonorCell parameter" << std::endl;
				}
			}
			else if (parameterName == "alpha")
			{
				alpha = atof(parameterValue.c_str());
			}
			else if (parameterName == "gamma")
			{
				gamma = atof(parameterValue.c_str());
			}
			else if (parameterName == "tau")
			{
				tau = atof(parameterValue.c_str());
			}
			else if (parameterName == "maximumDt")
			{
				maximumDt = atof(parameterValue.c_str());
			}
			else if (parameterName == "pressureSolver")
			{
				pressureSolver = parameterValue.c_str();
			}
			else if (parameterName == "omega")
			{
				omega = atof(parameterValue.c_str());
			}
			else if (parameterName == "epsilon")
			{
				epsilon = atof(parameterValue.c_str());
			}
			else if (parameterName == "maximumNumberOfIterations")
			{
				maximumNumberOfIterations = int(atof(parameterValue.c_str()));
			}
			else if (parameterName == "uInit")
			{
				uInit_ = atof(parameterValue.c_str());
			}
			else if (parameterName == "vInit")
			{
				vInit_ = atof(parameterValue.c_str());
			}
			else if (parameterName == "pInit")
			{
				pInit_ = atof(parameterValue.c_str());
			}
			else if (parameterName == "tInit")
			{
				TInit_ = atof(parameterValue.c_str());
			}
			else if (parameterName == "geometryFile")
			{
				geometryFile = parameterValue.c_str();
			}
			else if (parameterName == "outputFileEveryDt")
			{
				outputFileEveryDt = atof(parameterValue.c_str());
			}
			else if (parameterName == "participantName")
			{
				participantName= parameterValue.c_str();
			}
			else if (parameterName == "meshName")
			{
				meshName = parameterValue.c_str();
			}
			else if (parameterName == "preciceConfigFile")
			{
				preciceConfigFile = parameterValue.c_str();
			}
			else if (parameterName == "readDataName")
			{
				readDataName = parameterValue.c_str();
			}
			else if (parameterName == "writeDataName")
			{
				writeDataName = parameterValue.c_str();
			}
			else
			{
				std::cout << "Unknown parameter name: " << parameterName
						<< std::endl;
			}

		}

//    ...
//      else if (parameterName == "nCellsY")
//      {
//        nCells[1] = atoi(value.c_str());
//      }
//    ...

		// at the end of the file break for loop
		if (file.eof())
			break;

		// print line
		// std::cout << "line " << lineNo << ": " << line << std::endl;
	}

	loadGeometryFile();
};

void Settings::loadGeometryFile() {
	int nCellx;
	int nCelly;

	int TPD_count = 0;
	int TPN_count = 0;

	//geometry file
	// open file
	std::ifstream file(geometryFile.c_str(), std::ios::in);

	// check if file is open
	if (!file.is_open()) {
		std::cout << "Could not open geometry file " << geometryFile << std::endl;
		return;
	}

	// loop over lines of file
	for (int lineNo = 0;; lineNo++) {
		// read line
		std::string line;
		getline(file, line);
		if (line.find_first_of(" \t") != std::string::npos) {
			line.erase(0, line.find_first_not_of(" \t"));
		}

		std::string parameterName = line.substr(0, line.find_first_of('=')); //-1 ??
		// delete spaces after parameterName
		if (parameterName.find_first_of(" \t") != std::string::npos)
		{
			parameterName.erase(parameterName.find_first_of(" \t"));
		}
		std::string parameterValue = line.substr(line.find_first_of('=') + 1);

		// // delete spaces before parameterValue
		// if (parameterValue.find_first_of(" \t") != std::string::npos)
		// {
		// }
		// // delete commands
		// if (parameterValue.find_first_of("#") != std::string::npos)
		// {
		// 	parameterValue.erase(parameterValue.find_first_of("#"));
		// }

		// // delete spaces after parameterValue
		// if (parameterValue.find_first_of(" \t") != std::string::npos)
		// {
		// 	parameterValue.erase(parameterValue.find_first_of(" \t"));
		// }

		if (parameterName == "physicalSizeX")
		{
			physicalSize[0] = atof(parameterValue.c_str());
		}
		else if (parameterName == "physicalSizeY")
		{
			physicalSize[1] = atof(parameterValue.c_str());
		}
		else if (parameterName == "nCellsX")
		{
			nCells[0] = int(atof(parameterValue.c_str()));
		}
		else if (parameterName == "nCellsY")
		{
			nCells[1] = int(atof(parameterValue.c_str()));
			//initializing data structures
			nCellx = nCells[0] + 2;
			nCelly = nCells[1] + 2;
			std::array<int,2> nCellsGeometry = {nCellx, nCelly};
			geometryPVString_ = std::make_shared<Array2D>(nCellsGeometry); //nur Parameter, Rest egal
			geometryPVOrientation_ = std::make_shared<Array2D>(nCellsGeometry);
			geometryPV1_ = std::make_shared<Array2D>(nCellsGeometry);
			geometryPV2_ = std::make_shared<Array2D>(nCellsGeometry);
			geometryTString_ = std::make_shared<Array2D>(nCellsGeometry);
			geometryT1_ = std::make_shared<Array2D>(nCellsGeometry);
		}
		else if (parameterName == "xOrigin")
		{
			xOrigin = atof(parameterValue.c_str());
		}
		else if (parameterName == "yOrigin")
		{
			yOrigin = atof(parameterValue.c_str());
		}
		else if(parameterName == "Mesh")
		{
			for(int j = nCelly-1; j >= 0; j--) //(int i = nCelly-1; i >= 0; i--) // (int j = nCelly-1; j >= 0; j--)
			{
				getline(file, line);
				for(int i = 0; i < nCellx; i++) //(int j = nCellx-1; j >= 0; j--) // (int i = 0; i < nCellx; i++)
				{
					std::string cellAll = line.substr(0,line.find_first_of(","));
					line.erase(0,line.find_first_of(",")+1);
					if(cellAll.at(0) == 'F')
					{
						geometryPVString_->operator()(i,j) = -1;
						geometryTString_->operator()(i,j) = -1; // to avoid confilcts between TD and unassigned.
					}
					else
					{
						// determining the orientation of the boundary: (1,2,3,4 = left, right, lower, upper; 5,6,7,8 = lower-left, upper-left, lower-right, upper-right)
						// domain boundary cells
						if (i == 0) // && (j == 0 || geometryPVString_->operator()(i+1,j) == -1))
						{
							if (j == 0)
							{
								geometryPVOrientation_->operator()(i,j) = 8;
							}
							else if (j == nCelly-1)
							{
								geometryPVOrientation_->operator()(i,j) = 7;
							}
							else
							{
								geometryPVOrientation_->operator()(i,j) = 2;
							}
						}
						else if (i == nCellx-1) // && geometryPVString_->operator()(i-1,j) == -1)
						{
							if (j == 0)
							{
								geometryPVOrientation_->operator()(i,j) = 6;
							}
							else if (j == nCelly-1)
							{
								geometryPVOrientation_->operator()(i,j) = 5;
							}
							else
							{
								geometryPVOrientation_->operator()(i,j) = 1;
							}
						}
						else if (j == 0) // && geometryPVString_->operator()(i,j+1) == -1)
						{
							geometryPVOrientation_->operator()(i,j) = 4;
						}
						else if (j == nCelly-1) // && geometryPVString_->operator()(i,j-1) == -1)
						{
							geometryPVOrientation_->operator()(i,j) = 3;
						}
						else // inner cells
						{
							if (geometryPVString_->operator()(i-1,j) == -1) // left is fluid
							{
								geometryPVOrientation_->operator()(i,j) = 1;
							}
							if (geometryPVString_->operator()(i+1,j) == -1) // right
							{
								geometryPVOrientation_->operator()(i,j) = 2;
							}
							if (geometryPVString_->operator()(i,j-1) == -1) // lower
							{
								if (geometryPVOrientation_->operator()(i,j) == 1) // lower-left
								{
									geometryPVOrientation_->operator()(i,j) = 5;
								}
								else if (geometryPVOrientation_->operator()(i,j) == 2) // lower-right
								{
									geometryPVOrientation_->operator()(i,j) = 7;
								}
								else
								{
									geometryPVOrientation_->operator()(i,j) = 3;
								}
							}
							if (geometryPVString_->operator()(i,j+1) == -1) // upper
							{
								if (geometryPVOrientation_->operator()(i,j) == 1) // upper-left
								{
									geometryPVOrientation_->operator()(i,j) = 6;
								}
								else if (geometryPVOrientation_->operator()(i,j) == 2) // upper-right
								{
									geometryPVOrientation_->operator()(i,j) = 8;
								}
								else
								{
									geometryPVOrientation_->operator()(i,j) = 4;
								}
							}
						}

						std::string cellPressure = cellAll.substr(0,cellAll.find_first_of(";"));
						cellAll.erase(0, cellAll.find_first_of(";")+1);
						std::string cellTemperature = cellAll.substr(0);
						std::string cellPressureTyp = cellPressure.substr(0,cellPressure.find_first_of(":"));
						cellPressure.erase(0,cellPressure.find_first_of(":")+1);
						if(cellPressureTyp == "NSW")
						{
							geometryPVString_->operator()(i,j) = 0;
						}
						else if(cellPressureTyp == "SLW")
						{
							geometryPVString_->operator()(i,j) = 1;
						}
						else if(cellPressureTyp == "IN")
						{
							geometryPVString_->operator()(i,j) = 2;
							std::string cellPressure1 = cellPressure.substr(0,cellPressure.find_first_of(":"));
							cellPressure.erase(0,cellPressure.find_first_of(":")+1);
							std::string cellPressure2 = cellPressure.substr(0);
							geometryPV1_->operator()(i,j) = atof(cellPressure1.c_str());
							geometryPV2_->operator()(i,j) = atof(cellPressure2.c_str());
						}
						else if(cellPressureTyp == "OUT")
						{
							geometryPVString_->operator()(i,j) = 3;
						}
						else if(cellPressureTyp == "PR")
						{
							geometryPVString_->operator()(i,j) = 4;
							std::string cellPressure1 = cellPressure.substr(0,cellPressure.find_first_of(":"));
							cellPressure.erase(0,cellPressure.find_first_of(":")+1);
							geometryPV1_->operator()(i,j) = atof(cellPressure1.c_str());
						}
						else if (cellPressureTyp == "S")
						{
							geometryPVString_->operator()(i,j) = 5;
						}
						else
						{
							std::cout << "Unknow Cell Type!" << std::endl;
						}

						// New Temperature boundary conditions. Testing required.
						if (cellTemperature.find_first_of(",") <= cellTemperature.find_first_of(":"))
						{
							std::string cellTemperatureType = cellTemperature.substr(0,cellTemperature.find_first_of(","));
							if (cellTemperatureType == "TPD")
							{
								geometryTString_->operator()(i,j) = 2;
								TPD_count += 1;
								std::cout << "TPD count " << TPD_count << std::endl;

							}
							else if (cellTemperatureType == "TPN")
							{
								geometryTString_->operator()(i,j) = 3;
								TPN_count += 1;
								std::cout << "TPN count " << TPN_count << std::endl;

							}
							else
							{
								geometryTString_->operator()(i,j) = -1; // to avoid confilcts between TD and unassigned.
							}
						}
						else
						{
							std::string cellTemperatureType = cellTemperature.substr(0,cellTemperature.find_first_of(":"));
							cellTemperature.erase(0,cellTemperature.find_first_of(":")+1);
							if(cellTemperatureType == "TD")
							{
								geometryTString_->operator()(i,j) = 0; // waring: equal unassigned
								std::string cellTemperature1 = cellTemperature.substr(0);
								geometryT1_->operator()(i,j) = atof(cellTemperature1.c_str());
							}
							else if(cellTemperatureType == "TN")
							{
								geometryTString_->operator()(i,j) = 1;
								std::string cellTemperature1 = cellTemperature.substr(0);
								geometryT1_->operator()(i,j) = atof(cellTemperature1.c_str());
							}
							else
							{
								geometryTString_->operator()(i,j) = -1; // to avoid confilcts between TD and unassigned.
							}

						}
					}

				}
			}
		}
		// at the end of the file break for loop
		if (file.eof())
			break;

		// print line
		// std::cout << "line " << lineNo << ": " << line << std::endl;

	}

	// new preCICE related stuff: Getting interface (vertex) coordinates etc.
	if (TPD_count > 0 && TPN_count > 0)
	{
		std::cout << "This should not happen: TPD and PTN boundaries detected!" << std::endl;
	}

	//computing meshWidth
	double dx = physicalSize[0]/nCells[0];
	double dy = physicalSize[1]/nCells[1];

	vertexSize = std::max(TPD_count, TPN_count);
	std::cout << "vertexSize" << vertexSize << std::endl;

	vertex_i = std::vector<int> (vertexSize,0);
	vertex_j = std::vector<int> (vertexSize,0);
	std::vector<double> vertex_x(vertexSize,0);
	std::vector<double> vertex_y(vertexSize,0);
    // coords = new double[vertexSize*2];
	coords.resize(vertexSize*2);

	// std::vector<int> orientation_(vertexSize,0); // 0,1,2,3 = left,right,lower,upper
	orientation_.resize(vertexSize);
	int vertex_index = 0;
	for (int i = 0; i < nCells[0]+2; i++)
	{
		for (int j = 0; j < nCells[1]+2; j++)
		{
			if (geometryTString_->operator()(i,j) == 2 || geometryTString_->operator()(i,j) == 3)
			{
				vertex_i.at(vertex_index) = i;
				vertex_j.at(vertex_index) = j;
				vertex_x.at(vertex_index) = xOrigin + i*dx; // + .5*dx depending on the orientation of the boundary
				vertex_y.at(vertex_index) = yOrigin + j*dy; // + .5*dy depending on the orientation of the boundary
				vertex_index++ ;
			}
		}
	}
	for (int v = 0; v < vertexSize; v++)
	{
		orientation_[v] = geometryPVOrientation_->operator()(vertex_i.at(v), vertex_j.at(v));
		switch (orientation_[v])
		{
			case 0: vertex_x.at(v) -= .5*dx; break;
			case 1: vertex_x.at(v) += .5*dx; break;
			case 2: vertex_y.at(v) -= .5*dy; break;
			case 3:	vertex_y.at(v) += .5*dy; break;
			default: std::cout << "This shouldn't happen (unknown orientation): " << orientation_[v] << vertex_i.at(v) << vertex_j.at(v) <<std::endl; break;
		}
		coords.at(2*v) = vertex_x.at(v);
		coords.at(2*v+1) =vertex_y.at(v);

	}
}

void Settings::printSettings()
{
	std::cout << "Settings: " << std::endl
	<< "  physicalSize: " << physicalSize[0] << " x " << physicalSize[1] << ", nCells: " << nCells[0] << " x " << nCells[1] << std::endl
	<< "  endTime: " << endTime << " s, re: " << re << ", prandtl: " << prandtl << ", beta: " << beta << ", g: (" << g[0] << "," << g[1] << "), tau: "
	<< tau << ", maximum dt: " << maximumDt << std::endl
	<< "  dirichletBC: bottom: (" << dirichletBcBottom[0] << "," << dirichletBcBottom[1] << ")"
	<< ", top: (" << dirichletBcTop[0] << "," << dirichletBcTop[1] << ")"
	<< ", left: (" << dirichletBcLeft[0] << "," << dirichletBcLeft[1] << ")"
	<< ", right: (" << dirichletBcRight[0] << "," << dirichletBcRight[1] << ")" << std::endl
	<< "  useDonorCell: " << std::boolalpha << useDonorCell << ", alpha: " << alpha << std::endl
	<< "  pressureSolver: " << pressureSolver << ", omega: " << omega << ", epsilon: " << epsilon << ", maximumNumberOfIterations: " << maximumNumberOfIterations << std::endl;
}
