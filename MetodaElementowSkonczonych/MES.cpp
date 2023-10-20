#include "MES.h"

void GlobalData::showGlobalData()
{
	std::cout << SimulationTime << std::endl;
	std::cout << SimulationStepTime << std::endl;
	std::cout << Conductivity << std::endl;
	std::cout << Alfa << std::endl;
	std::cout << Tot << std::endl;
	std::cout << InitialTemp << std::endl;
	std::cout << Density << std::endl;
	std::cout << SpecificHeat << std::endl;
	std::cout << NodesNumber << std::endl;
	std::cout << ElementsNumber << std::endl;
}

void ElementUniwersalny::init(int n)
{
	int nPktCalk = pow(n, 2);

	matdEta = new double* [nPktCalk];
	matdKsi = new double* [nPktCalk];

	for (int i = 0; i < nPktCalk; i++)
	{
		matdEta[i] = new double[4];
		matdKsi[i] = new double[4];
	}


	for (int i = 0; i < nPktCalk; i++)
	{
		for (int j = 0; j < 4; j++)
		{

		}
	}
}

GlobalData* readMesh(std::string fileSrc)
{
	GlobalData* output = new GlobalData;
	std::string line;
	std::string value;

	std::ifstream file;
	file.open(fileSrc);
	if (file.good())
	{
		getline(file, line);
		value = line.substr(line.find_last_of(' ') + 1, line.length());
		output->SimulationTime = stoi(value);

		getline(file, line);
		value = line.substr(line.find_last_of(' ') + 1, line.length());
		output->SimulationStepTime = stoi(value);

		getline(file, line);
		value = line.substr(line.find_last_of(' ') + 1, line.length());
		output->Conductivity = stoi(value);

		getline(file, line);
		value = line.substr(line.find_last_of(' ') + 1, line.length());
		output->Alfa = stoi(value);

		getline(file, line);
		value = line.substr(line.find_last_of(' ') + 1, line.length());
		output->Tot = stoi(value);

		getline(file, line);
		value = line.substr(line.find_last_of(' ') + 1, line.length());
		output->InitialTemp = stoi(value);

		getline(file, line);
		value = line.substr(line.find_last_of(' ') + 1, line.length());
		output->Density = stoi(value);

		getline(file, line);
		value = line.substr(line.find_last_of(' ') + 1, line.length());
		output->SpecificHeat = stoi(value);

		getline(file, line);
		value = line.substr(line.find_last_of(' ') + 1, line.length());
		output->NodesNumber = stoi(value);

		getline(file, line);
		value = line.substr(line.find_last_of(' ') + 1, line.length());
		output->ElementsNumber = stoi(value);
	}

	return output;
}
