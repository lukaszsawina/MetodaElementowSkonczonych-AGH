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
	nPkt = n;
	int nn = n;
	double* xn = new double[nn];

	switch (n)
	{
	case 2:
		xn[0] = -0.577350;
		xn[1] = 0.577350;
		break;
	case 3:
		xn[0] = -0.774597;
		xn[1] = 0;
		xn[2] = 0.774597;
		break;
	case 4:
		xn[0] = -0.861136;
		xn[1] = -0.339981;
		xn[2] = 0.339981;
		xn[3] = 0.861136;
		break;
	case 5:
		xn[0] = -0.906180;
		xn[1] = -0.538469;
		xn[2] = 0;
		xn[3] = 0.538469;
		xn[4] = 0.906180;
		break;
	}

	int nPktCalk = pow(n, 2);

	matdEta = new double* [nPktCalk];
	matdKsi = new double* [nPktCalk];

	for (int i = 0; i < nPktCalk; i++)
	{
		matdEta[i] = new double[4];
		matdKsi[i] = new double[4];
	}

	for (int i = 0, k =0; i < n; i++)
	{
		for (int j = 0; j < n; j++, k++)
		{
			matdKsi[k][0] = -0.25 * (1 - xn[i]);
			matdKsi[k][1] = 0.25 * (1 - xn[i]);
			matdKsi[k][2] = 0.25 * (1 + xn[i]);
			matdKsi[k][3] = -0.25 * (1 + xn[i]);

			matdEta[k][0] = -0.25 * (1 - xn[j]);
			matdEta[k][1] = -0.25 * (1 + xn[j]);
			matdEta[k][2] = 0.25 * (1 + xn[j]);
			matdEta[k][3] = 0.25 * (1 - xn[j]);
		}
	}
	delete[] xn;
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
