#include "Mesh.h"

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
		std::cout << "ID: " << nodes[i].ID << " X: " << nodes[i].x << " Y: " << nodes[i].y << " BC: " << nodes[i].BC << std::endl;	
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

void Mesh::calcHForElements(const ElementUniwersalny& elUni)
{
	for (int i = 0; i < globalData->ElementsNumber; i++)
	{
		std::cout << std::endl << "Element: " << i+1 << std::endl;
		double* x = new double[4];
		double* y = new double[4];

		for (int j = 0; j < 4; j++)
		{
			int n = elements[i].ID_wezlow[j];
			x[j] = nodes[n-1].x;
			y[j] = nodes[n-1].y;

			std::cout << "PC " << i+1 << "(x: " << x[j] << " ; " << y[j] << ")" << std::endl;
		}

		elements[i].calcH(x, y, elUni);

		delete[] x;
		delete[] y;
	}
}

void Mesh::calcHBCForElements(const ElementUniwersalny& elUni)
{
	for (int i = 0; i < globalData->ElementsNumber; i++)
	{
		std::cout << std::endl <<  "Element: " << i + 1 << std::endl;
		double* x = new double[4];
		double* y = new double[4];
		int* bc = new int[4];

		for (int j = 0; j < 4; j++)
		{
			int n = elements[i].ID_wezlow[j];
			x[j] = nodes[n - 1].x;
			y[j] = nodes[n - 1].y;	
			bc[j] = nodes[n - 1].BC;
		}

		elements[i].calcHBC(x, y, bc, elUni);

		delete[] x;
		delete[] y;
		delete[] bc;
	}
}

void Mesh::calcVectorPForElements(const ElementUniwersalny& elUni)
{
	for (int i = 0; i < globalData->ElementsNumber; i++)
	{
		std::cout << std::endl << "Element: " << i + 1 << std::endl;
		double* x = new double[4];
		double* y = new double[4];
		int* bc = new int[4];

		for (int j = 0; j < 4; j++)
		{
			int n = elements[i].ID_wezlow[j];
			x[j] = nodes[n - 1].x;
			y[j] = nodes[n - 1].y;
			bc[j] = nodes[n - 1].BC;
		}

		elements[i].calcVectorP(x, y, bc, elUni);

		delete[] x;
		delete[] y;
		delete[] bc;
	}
}

void Mesh::calcCForElements(const ElementUniwersalny& elUni)
{
	for (int i = 0; i < globalData->ElementsNumber; i++)
	{
		std::cout << std::endl << "Element: " << i + 1 << std::endl;
		double* x = new double[4];
		double* y = new double[4];

		for (int j = 0; j < 4; j++)
		{
			int n = elements[globalData->ElementsNumber - 1 - i].ID_wezlow[j];
			x[j] = nodes[n - 1].x;
			y[j] = nodes[n - 1].y;
		}

		elements[i].calcC(x, y, elUni);

		delete[] x;
		delete[] y;
	}
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

		getline(file, line);
		while (line != "*BC")
			getline(file, line);

		std::vector<std::string> elements;
		std::string element;

		getline(file, line);

		std::istringstream BCstream(line);

		while (std::getline(BCstream, element, ','))
		{
			size_t start = element.find_first_not_of(" ");
			size_t end = element.find_last_not_of(" ");
			elements.push_back(element.substr(start, end - start + 1));
		}

		for (int i = 0; i < elements.size(); i++)
		{
			for (int j = 0; j < globalData->NodesNumber; j++)
				if (output[j].ID == std::stoi(elements[i]))
					output[j].BC = 1;

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

Mesh::~Mesh()
{
	for (int i = 0; i < globalData->ElementsNumber; i++)
	{
		if (elements[i].H != nullptr)
		{
			for (int j = 0; j < 4; j++)
				delete[] elements[i].H[j];
			delete[] elements[i].H;
		}

		if (elements[i].HBC != nullptr)
		{
			for (int j = 0; j < 4; j++)
				delete[] elements[i].HBC[j];
			delete[] elements[i].HBC;
		}

		if(elements[i].VectorP != nullptr)
			delete[] elements[i].VectorP;
	}

	delete[] elements;
	delete[] nodes;

	delete globalData;
}

// Funkcja znajduj¹ca wartoœæ minimaln¹ i maksymaln¹ w tablicy
std::pair<double, double> znajdz_min_i_max(const double* tablica, int rozmiar) {
	if (rozmiar == 0 || tablica == nullptr) {
		// Zwracamy wartoœci domyœlne dla pustej tablicy lub nullptr
		return std::make_pair(std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity());
	}

	double min_wartosc = std::numeric_limits<double>::infinity();
	double max_wartosc = -std::numeric_limits<double>::infinity();

	for (int i = 0; i < rozmiar; ++i) {
		if (tablica[i] < min_wartosc) {
			min_wartosc = tablica[i];
		}
		if (tablica[i] > max_wartosc) {
			max_wartosc = tablica[i];
		}
	}

	return std::make_pair(min_wartosc, max_wartosc);
}

double* calcTemperatureForStep(Mesh& mesh, double* tempV)
{
	int nNodes = mesh.globalData->NodesNumber;

	double** agregatH = new double* [nNodes];

	for (int i = 0; i < nNodes; i++)
	{
		agregatH[i] = new double[nNodes];

		for (int j = 0; j < nNodes; j++)
			agregatH[i][j] = 0;
	}

	double** agregatC = new double* [nNodes];

	for (int i = 0; i < nNodes; i++)
	{
		agregatC[i] = new double[nNodes];

		for (int j = 0; j < nNodes; j++)
			agregatC[i][j] = 0;
	}

	double* agregatP = new double[nNodes];

	for (int i = 0; i < nNodes; i++)
	{
		agregatP[i] = 0;
	}

	//Liczenie agregatu H
	int nElements = mesh.globalData->ElementsNumber;

	for (int e = 0; e < nElements; e++)
	{
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				agregatH[mesh.elements[e].ID_wezlow[i] - 1][mesh.elements[e].ID_wezlow[j] - 1] += mesh.elements[e].H[i][j];
				agregatH[mesh.elements[e].ID_wezlow[i] - 1][mesh.elements[e].ID_wezlow[j] - 1] += mesh.elements[e].HBC[i][j];
				agregatH[mesh.elements[e].ID_wezlow[i] - 1][mesh.elements[e].ID_wezlow[j] - 1] += mesh.elements[e].C[i][j] / mesh.globalData->SimulationStepTime;
				agregatC[mesh.elements[e].ID_wezlow[i] - 1][mesh.elements[e].ID_wezlow[j] - 1] += mesh.elements[e].C[i][j] / mesh.globalData->SimulationStepTime;
			}
		}
	}

	//std::cout << std::endl << "Agregat H" << std::endl;

	//for (int i = 0; i < nNodes; i++)
	//{
	//	for (int j = 0; j < nNodes; j++)
	//		std::cout << agregatC[i][j] << "\t";
	//	std::cout << std::endl;
	//}

	double* agregatBC = new double[nNodes];

	for (int i = 0; i < nNodes; i++)
	{
		agregatBC[i] = 0;
	}

	for (int i = 0; i < nNodes; i++)
	{
		double temp = 0;
		for (int j = 0; j < nNodes; j++)
		{
			temp += agregatC[i][j] * tempV[j];
		}
		agregatBC[i] = temp;
	}

	//Liczenie agregatu P
	for (int e = 0; e < nElements; e++)
	{
		for (int i = 0; i < 4; i++)
		{
			agregatP[mesh.elements[e].ID_wezlow[i] - 1] += mesh.elements[e].VectorP[i];
		}
	}

	for (int i = 0; i < nNodes; i++)
	{
		agregatBC[i] = (agregatBC[i] + agregatP[i]);
	}

	//std::cout << std::endl << "Agregat BC" << std::endl;

	//for (int i = 0; i < nNodes; i++)
	//{
	//	std::cout << agregatBC[i] << std::endl;
	//}

	double* output = ElimGauss(agregatH, agregatBC, nNodes);

	return output;
}

double* Mesh::calcTemperature(Mesh& mesh, const ElementUniwersalny& elUni)
{
	mesh.calcHForElements(elUni);
	mesh.calcHBCForElements(elUni);
	mesh.calcVectorPForElements(elUni);
	mesh.calcCForElements(elUni);

	double* tempV = new double[mesh.globalData->NodesNumber];

	for (int i = 0; i < mesh.globalData->NodesNumber; i++)
	{
		tempV[i] = mesh.globalData->InitialTemp;
	}

	int nstep = mesh.globalData->SimulationTime / mesh.globalData->SimulationStepTime;

	for (int i = 0; i < nstep; i++)
	{
		tempV = calcTemperatureForStep(mesh, tempV);

		std::cout << std::endl;
		std::cout << "Time: " << mesh.globalData->SimulationStepTime * (i + 1) << std::endl;

		std::pair<double, double> wyniki = znajdz_min_i_max(tempV, mesh.globalData->NodesNumber);
		std::cout << "Min: " << wyniki.first << "\t" << wyniki.second;
		/*for (int j = 0; j < mesh.globalData->NodesNumber; j++)
		{
			std::cout << tempV[j] << " ";
		}*/
		std::cout << std::endl;
	}

	return tempV;
}

