#include "MES.h"

void Mesh::showGlobalData()
{
	std::cout << globalData->SimulationTime << std::endl;
	std::cout << globalData->SimulationStepTime << std::endl;
	std::cout << globalData->Conductivity << std::endl;
	std::cout << globalData->Alfa << std::endl;
	std::cout << globalData->Tot << std::endl;
	std::cout << globalData->InitialTemp << std::endl;
	std::cout << globalData->Density << std::endl;
	std::cout << globalData->SpecificHeat << std::endl;
	std::cout << globalData->NodesNumber << std::endl;
	std::cout << globalData->ElementsNumber << std::endl;
}

void Mesh::showNodes()
{
	std::cout << std::setprecision(10);
	for (int i = 0; i < globalData->NodesNumber; i++)
		std::cout << "ID: " << nodes[i].ID << " X: " << nodes[i].x << " Y: " << nodes[i].y << std::endl;
		
}

void Mesh::showElements()
{
	std::cout << std::setprecision(10);
	for (int i = 0; i < globalData->ElementsNumber; i++)
	{
		std::cout << "Element ID: " << elements[i].ID << std::endl;
		for (int j = 0; j < 4; j++)
			std::cout << elements[i].ID_wezlow[j] << "\t";
		std::cout << std::endl;
	}
}

void ElementUniwersalny::init(int n)
{
	nPkt = n;
	double* G_X = initXGauss(n);

	int nPktCalk = pow(n, 2);

	matdEta = new double* [nPktCalk];
	matdKsi = new double* [nPktCalk];

	for (int i = 0; i < nPktCalk; i++)
	{
		matdEta[i] = new double[4];
		matdKsi[i] = new double[4];
	}

	//Liczenie wartoœci macierzy ksi i eta
	for (int i = 0, k =0; i < n; i++)
	{
		for (int j = 0; j < n; j++, k++)
		{
			matdKsi[k][0] = -0.25 * (1 - G_X[i]);
			matdKsi[k][1] = 0.25 * (1 - G_X[i]);
			matdKsi[k][2] = 0.25 * (1 + G_X[i]);
			matdKsi[k][3] = -0.25 * (1 + G_X[i]);

			matdEta[k][0] = -0.25 * (1 - G_X[j]);
			matdEta[k][1] = -0.25 * (1 + G_X[j]);
			matdEta[k][2] = 0.25 * (1 + G_X[j]);
			matdEta[k][3] = 0.25 * (1 - G_X[j]);
		}
	}
	delete[] G_X;
}

double** ElementUniwersalny::H(double* x, double* y)
{
	std::vector<double**> macierzeJakobiegoPunktow;
	std::vector<double> detJPunktow;

	for (int i = 0; i < pow(nPkt, 2); i++)
	{
		double** macierzJakobiego = new double* [2];

		for (int i = 0; i < 2; i++)
		{
			macierzJakobiego[i] = new double[2];
		}

		//Liczenie macierzy jakobiego dla pc[i]
		macierzJakobiego[0][0] = matdEta[i][0] * y[0] + matdEta[i][1] * y[1] + matdEta[i][2] * y[2] + matdEta[i][3] * y[3];
		macierzJakobiego[0][1] = -1 * (matdKsi[i][0] * y[0] + matdKsi[i][1] * y[1] + matdKsi[i][2] * y[2] + matdKsi[i][3] * y[3]);
		macierzJakobiego[1][0] = -1 * (matdEta[i][0] * x[0] + matdEta[i][1] * x[1] + matdEta[i][2] * x[2] + matdEta[i][3] * x[3]);
		macierzJakobiego[1][1] = matdKsi[i][0] * x[0] + matdKsi[i][1] * x[1] + matdKsi[i][2] * x[2] + matdKsi[i][3] * x[3];

		macierzeJakobiegoPunktow.push_back(macierzJakobiego);

		double detJ = macierzJakobiego[0][0] * macierzJakobiego[1][1] - (macierzJakobiego[1][0] * macierzJakobiego[0][1]);
		detJPunktow.push_back(detJ);
	}

	std::cout << "Macierze jakobiego punktów" << std::endl;
	for (int i = 0; i < macierzeJakobiegoPunktow.size(); i++)
	{
		std::cout << "PC" << i + 1 << std::endl;

		for (int j = 0; j < 2; j++)
		{
			for (int k = 0; k < 2; k++)
				std::cout << macierzeJakobiegoPunktow[i][j][k] << "\t";
			std::cout << std::endl;
		}
		std::cout << std::endl;
	}

	//Przemno¿enie przez 1/detJ
	for (int i = 0; i < macierzeJakobiegoPunktow.size(); i++)
	{
		for (int j = 0; j < 2; j++)
			for (int k = 0; k < 2; k++)
				macierzeJakobiegoPunktow[i][j][k] *= (1 / detJPunktow[i]);
	}

	std::cout << "Macierze jakobiego punktów przemno¿ona przez 1/detJ" << std::endl;
	for (int i = 0; i < macierzeJakobiegoPunktow.size(); i++)
	{
		std::cout << "PC" << i + 1 << std::endl;

		for (int j = 0; j < 2; j++)
		{
			for (int k = 0; k < 2; k++)
				std::cout << macierzeJakobiegoPunktow[i][j][k] << "\t";
			std::cout << std::endl;
		}
		std::cout << std::endl;
	}

	double** matdx;
	double** matdy;

	matdx = new double* [4];
	matdy = new double* [4];

	for (int i = 0; i < 4; i++)
	{
		matdx[i] = new double[4];
		matdy[i] = new double[4];
	}

	//Liczenie wartoœci macierzy dx i dy
	for (int i = 0; i < 4; i++)
	{
		matdx[i][0] = macierzeJakobiegoPunktow[i][0][0] * matdKsi[i][0] + macierzeJakobiegoPunktow[i][0][1] * matdEta[i][0];
		matdx[i][1] = macierzeJakobiegoPunktow[i][0][0] * matdKsi[i][1] + macierzeJakobiegoPunktow[i][0][1] * matdEta[i][1];
		matdx[i][2] = macierzeJakobiegoPunktow[i][0][0] * matdKsi[i][2] + macierzeJakobiegoPunktow[i][0][1] * matdEta[i][2];
		matdx[i][3] = macierzeJakobiegoPunktow[i][0][0] * matdKsi[i][3] + macierzeJakobiegoPunktow[i][0][1] * matdEta[i][3];

		matdy[i][0] = macierzeJakobiegoPunktow[i][1][0] * matdKsi[i][0] + macierzeJakobiegoPunktow[i][1][1] * matdEta[i][0];
		matdy[i][1] = macierzeJakobiegoPunktow[i][1][0] * matdKsi[i][1] + macierzeJakobiegoPunktow[i][1][1] * matdEta[i][1];
		matdy[i][2] = macierzeJakobiegoPunktow[i][1][0] * matdKsi[i][2] + macierzeJakobiegoPunktow[i][1][1] * matdEta[i][2];
		matdy[i][3] = macierzeJakobiegoPunktow[i][1][0] * matdKsi[i][3] + macierzeJakobiegoPunktow[i][1][1] * matdEta[i][3];
	}

	//Zwolnienei pamiêci z macierzyJakobiegoPunktów
	for (int i = 0; i < macierzeJakobiegoPunktow.size(); i++)
	{
		for (int j = 0; j < 2; j++)
			delete[] macierzeJakobiegoPunktow[i][j];

		delete[] macierzeJakobiegoPunktow[i];
	}

	std::cout << "matdx" << std::endl;
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
			std::cout << matdx[i][j] << " ";
		std::cout << std::endl;
	}

	std::cout << "matdy" << std::endl;
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
			std::cout << matdy[i][j] << " ";
		std::cout << std::endl;
	}

	std::vector<double**> macierzeHPunktow;

	for (int p = 0; p < pow(nPkt, 2); p++)
	{
		double** Hpkt = new double* [4];

		for (int i = 0; i < 4; i++)
			Hpkt[i] = new double[4];

		//Liczenie macierzy H dla pc[i]
		for (int i = 0; i < 4; i++)
		{
			//Wartoœæ 30 nie jest sta³a!!! pewnie bêdzie zmieniana w dalszych programach
			Hpkt[i][0] = (matdx[p][i] * matdx[p][0] + matdy[p][i] * matdy[p][0]) * detJPunktow[p] * 30;
			Hpkt[i][1] = (matdx[p][i] * matdx[p][1] + matdy[p][i] * matdy[p][1]) * detJPunktow[p] * 30;
			Hpkt[i][2] = (matdx[p][i] * matdx[p][2] + matdy[p][i] * matdy[p][2]) * detJPunktow[p] * 30;
			Hpkt[i][3] = (matdx[p][i] * matdx[p][3] + matdy[p][i] * matdy[p][3]) * detJPunktow[p] * 30;
		}

		macierzeHPunktow.push_back(Hpkt);
	}

	std::cout << std::endl;

	std::cout << "Macierze H punktów" << std::endl;

	for (int i = 0; i < macierzeHPunktow.size(); i++)
	{
		std::cout << "PC" << i + 1 << std::endl;

		for (int j = 0; j < 4; j++)
		{
			for (int k = 0; k < 4; k++)
				std::cout << macierzeHPunktow[i][j][k] << "\t";
			std::cout << std::endl;
		}
		std::cout << std::endl;
	}


	double** H = new double* [4];
	for (int i = 0; i < 4; i++)
		H[i] = new double[4];

	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
			H[i][j] = 0;
	}

	double* G_W = initWGauss(nPkt);

	//Liczenie macierzy H dla elementu
	for (int i = 0, p = 0; i < nPkt; i++)
	{
		for (int j = 0; j < nPkt; j++, p++)
		{
			for (int n = 0; n < 4; n++)
				for (int m = 0; m < 4; m++)
				{
					macierzeHPunktow[p][n][m] *= (G_W[i] * G_W[j]);
					H[n][m] += macierzeHPunktow[p][n][m];
				}
		}
	}

	//Zwolnienie pamiêci dla macierzy H punktów
	delete[] G_W;
	for (int i = 0; i < macierzeHPunktow.size(); i++)
	{
		for (int j = 0; j < 4; j++)
			delete[] macierzeHPunktow[i][j];

		delete[] macierzeHPunktow[i];
	}

	std::cout << "Macierz H" << std::endl;
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
			std::cout << H[i][j] << "\t";
		std::cout << std::endl;
	}

	return H;
}

GlobalData* Mesh::readMeshGlobalData(std::string fileSrc)
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

Node* Mesh::readMeshNodes(std::string fileSrc)
{
	Node* output = new Node[globalData->NodesNumber];

	std::string line;

	std::ifstream file;
	file.open(fileSrc);
	if (file.good())
	{
		//Szukamy miejsca, gdzie zaczynaj¹ siê nody
		getline(file, line);
		while(line != "*Node")
			getline(file, line);

		for (int i = 0; i < globalData->NodesNumber; i++)
		{
			getline(file, line);
			std::string t = line;

			std::istringstream stream(line);
			std::vector<std::string> elements;

			std::string element;
			while (std::getline(stream, element, ','))
			{
				size_t start = element.find_first_not_of(" ");
				size_t end = element.find_last_not_of(" ");
				elements.push_back(element.substr(start, end - start + 1));
			}

			output[i].ID = stoi(elements[0]);
			output[i].x = stod(elements[1]);
			output[i].y = stod(elements[2]);
		}
	}

	return output;
}

Element* Mesh::readMeshElements(std::string fileSrc)
{
	Element* output = new Element[globalData->ElementsNumber];

	std::string line;

	std::ifstream file;
	file.open(fileSrc);
	if (file.good())
	{
		//Szukamy miejsca, gdzie zaczynaj¹ siê elementy
		getline(file, line);
		while (line.find("*Element"))
			getline(file, line);

		for (int i = 0; i < globalData->ElementsNumber; i++)
		{
			getline(file, line);
			std::string t = line;

			std::istringstream stream(line);
			std::vector<std::string> elements;

			std::string element;
			while (std::getline(stream, element, ','))
			{
				size_t start = element.find_first_not_of(" ");
				size_t end = element.find_last_not_of(" ");
				elements.push_back(element.substr(start, end - start + 1));
			}

			output[i].ID = stoi(elements[0]);

			for (int j = 0; j < 4; j++)
				output[i].ID_wezlow[j] = stod(elements[j+1]);
		}
	}

	return output;
}

void Mesh::readMeshFile(std::string fileSrc)
{
	globalData = readMeshGlobalData(fileSrc);
	nodes = readMeshNodes(fileSrc);
	elements = readMeshElements(fileSrc);
}