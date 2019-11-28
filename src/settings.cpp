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
//    while(line[0] == ' ' || line[0]== '\t')
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
				beta = int(atof(parameterValue.c_str()));
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
					std::cout << "nix funktioniert" << std::endl;
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
				uInit = atof(parameterValue.c_str());
			}
			else if (parameterName == "vInit")
			{
				vInit = atof(parameterValue.c_str());
			}
			else if (parameterName == "pInit")
			{
				pInit = atof(parameterValue.c_str());
			}
			else if (parameterName == "tInit")
			{
				tInit = atof(parameterValue.c_str());
			}
			else if (parameterName == "geometryFile")
			{
				geometryFile = parameterValue.c_str();
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
	//initializing data structures
	int nCellx = nCells[0] + 2;
	int nCelly = nCells[1] + 2;
	std::array<int,2> nCellsGeometry = {nCellx, nCelly};
	geometryPVString = std::make_shared<Array2D>(nCellsGeometry); //nur Parameter, Rest egal
	geometryPV1 = std::make_shared<Array2D>(nCellsGeometry);
	geometryPV2 = std::make_shared<Array2D>(nCellsGeometry);
	geometryTString = std::make_shared<Array2D>(nCellsGeometry);
	geometryT1 = std::make_shared<Array2D>(nCellsGeometry);


	//geometry file
	// open file
	std::ifstream file(geometryFile.c_str(), std::ios::in);
	std::cout <<  geometryFile << std::endl;

	// check if file is open
	if (!file.is_open()) {
		std::cout << "Could not  open geometry file " << geometryFile << "Komische Bugs." << std::endl;
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
		// delete spaces before parameterValue
		if (parameterValue.find_first_of(" \t") != std::string::npos)
		{
		parameterValue.erase(0,parameterValue.find_first_not_of(" \t"));
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
		else if (parameterName == "nCellsX")
		{
			nCells[0] = int(atof(parameterValue.c_str()));
		}
		else if (parameterName == "nCellsY")
		{
			nCells[1] = int(atof(parameterValue.c_str()));
		}
		 else if(parameterName == "Mesh")
		{
			for(int i = nCelly-1; i >= 0; i--)
			{
				getline(file, line);
				for(int j = nCellx-1; j >= 0; j--)
				{
					std::string cellAll = line.substr(0,line.find_first_of(","));
					line.erase(0,line.find_first_of(",")+1);
					if(cellAll.at(0) != 'F')
					{
						std::string cellPressure = cellAll.substr(0,cellAll.find_first_of(";"));
						cellAll.erase(0, cellAll.find_first_of(";")+1);
						std::string cellTemperature = cellAll.substr(0);
						std::string cellPressureTyp = cellPressure.substr(0,cellPressure.find_first_of(":"));
						cellPressure.erase(0,cellPressure.find_first_of(":")+1);
						if(cellPressureTyp == "NSW")
						{
							geometryPVString->operator()(i,j) = 0;
						}
						else if(cellPressureTyp == "SLW")
						{
							geometryPVString->operator()(i,j) = 1;
						}
						else if(cellPressureTyp == "IN")
						{
							geometryPVString->operator()(i,j) = 2;
							std::string cellPressure1 = cellPressure.substr(0,cellPressure.find_first_of(":"));
							cellPressure.erase(0,cellPressure.find_first_of(":")+1);
							std::string cellPressure2 = cellPressure.substr(0);
							geometryPV1->operator()(i,j) = atof(cellPressure1.c_str());
							geometryPV2->operator()(i,j) = atof(cellPressure2.c_str());
						}
						else if(cellPressureTyp == "OUT")
						{
							geometryPVString->operator()(i,j) = 3;
						}
						else if(cellPressureTyp == "PR")
						{
							geometryPVString->operator()(i,j) = 4;
							std::string cellPressure1 = cellPressure.substr(0,cellPressure.find_first_of(":"));
							cellPressure.erase(0,cellPressure.find_first_of(":")+1);
							geometryPV1->operator()(i,j) = atof(cellPressure1.c_str());
						}
						std::string cellTemperatureType = cellTemperature.substr(0,cellTemperature.find_first_of(":"));
						cellTemperature.erase(0,cellTemperature.find_first_of(":")+1);
						if(cellTemperatureType == "TD")
						{
							geometryTString->operator()(i,j) = 0;
							std::string cellTemperature1 = cellTemperature.substr(0);
							geometryT1->operator()(i,j) = atof(cellTemperature1.c_str());
						} 
						else if(cellTemperatureType == "TN")
						{
							geometryTString->operator()(i,j) = 1;
							std::string cellTemperature1 = cellTemperature.substr(0);
							geometryT1->operator()(i,j) = atof(cellTemperature1.c_str());
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
