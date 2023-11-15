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

void Element::calcH(double* x, double* y, ElementUniwersalny elUni)
{
	std::vector<double**> macierzeJakobiegoPunktow;
	std::vector<double> detJPunktow;

	for (int i = 0; i < pow(elUni.nPkt, 2); i++)
	{
		double** macierzJakobiego = new double* [2];

		for (int i = 0; i < 2; i++)
		{
			macierzJakobiego[i] = new double[2];
		}

		//Liczenie macierzy jakobiego dla pc[i]
		macierzJakobiego[0][0] = elUni.matdEta[i][0] * y[0] + elUni.matdEta[i][1] * y[1] + elUni.matdEta[i][2] * y[2] + elUni.matdEta[i][3] * y[3];
		macierzJakobiego[0][1] = -1 * (elUni.matdKsi[i][0] * y[0] + elUni.matdKsi[i][1] * y[1] + elUni.matdKsi[i][2] * y[2] + elUni.matdKsi[i][3] * y[3]);
		macierzJakobiego[1][0] = -1 * (elUni.matdEta[i][0] * x[0] + elUni.matdEta[i][1] * x[1] + elUni.matdEta[i][2] * x[2] + elUni.matdEta[i][3] * x[3]);
		macierzJakobiego[1][1] = elUni.matdKsi[i][0] * x[0] + elUni.matdKsi[i][1] * x[1] + elUni.matdKsi[i][2] * x[2] + elUni.matdKsi[i][3] * x[3];

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
		matdx[i][0] = macierzeJakobiegoPunktow[i][0][0] * elUni.matdKsi[i][0] + macierzeJakobiegoPunktow[i][0][1] * elUni.matdEta[i][0];
		matdx[i][1] = macierzeJakobiegoPunktow[i][0][0] * elUni.matdKsi[i][1] + macierzeJakobiegoPunktow[i][0][1] * elUni.matdEta[i][1];
		matdx[i][2] = macierzeJakobiegoPunktow[i][0][0] * elUni.matdKsi[i][2] + macierzeJakobiegoPunktow[i][0][1] * elUni.matdEta[i][2];
		matdx[i][3] = macierzeJakobiegoPunktow[i][0][0] * elUni.matdKsi[i][3] + macierzeJakobiegoPunktow[i][0][1] * elUni.matdEta[i][3];

		matdy[i][0] = macierzeJakobiegoPunktow[i][1][0] * elUni.matdKsi[i][0] + macierzeJakobiegoPunktow[i][1][1] * elUni.matdEta[i][0];
		matdy[i][1] = macierzeJakobiegoPunktow[i][1][0] * elUni.matdKsi[i][1] + macierzeJakobiegoPunktow[i][1][1] * elUni.matdEta[i][1];
		matdy[i][2] = macierzeJakobiegoPunktow[i][1][0] * elUni.matdKsi[i][2] + macierzeJakobiegoPunktow[i][1][1] * elUni.matdEta[i][2];
		matdy[i][3] = macierzeJakobiegoPunktow[i][1][0] * elUni.matdKsi[i][3] + macierzeJakobiegoPunktow[i][1][1] * elUni.matdEta[i][3];
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

	for (int p = 0; p < pow(elUni.nPkt, 2); p++)
	{
		double** Hpkt = new double* [4];

		for (int i = 0; i < 4; i++)
			Hpkt[i] = new double[4];

		//Liczenie macierzy H dla pc[i]
		for (int i = 0; i < 4; i++)
		{
			//Wartoœæ 30 nie jest sta³a!!! pewnie bêdzie zmieniana w dalszych programach
			Hpkt[i][0] = (matdx[p][i] * matdx[p][0] + matdy[p][i] * matdy[p][0]) * detJPunktow[p] * 25;
			Hpkt[i][1] = (matdx[p][i] * matdx[p][1] + matdy[p][i] * matdy[p][1]) * detJPunktow[p] * 25;
			Hpkt[i][2] = (matdx[p][i] * matdx[p][2] + matdy[p][i] * matdy[p][2]) * detJPunktow[p] * 25;
			Hpkt[i][3] = (matdx[p][i] * matdx[p][3] + matdy[p][i] * matdy[p][3]) * detJPunktow[p] * 25;
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


	double** outputH = new double* [4];
	for (int i = 0; i < 4; i++)
		outputH[i] = new double[4];

	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
			outputH[i][j] = 0;
	}

	double* G_W = initWGauss(elUni.nPkt);

	//Liczenie macierzy H dla elementu
	for (int i = 0, p = 0; i < elUni.nPkt; i++)
	{
		for (int j = 0; j < elUni.nPkt; j++, p++)
		{
			for (int n = 0; n < 4; n++)
				for (int m = 0; m < 4; m++)
				{
					macierzeHPunktow[p][n][m] *= (G_W[i] * G_W[j]);
					outputH[n][m] += macierzeHPunktow[p][n][m];
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
			std::cout << outputH[i][j] << "\t";
		std::cout << std::endl;
	}

	H = outputH;
}

void Element::calcHBC(double* x, double* y, double* BC, ElementUniwersalny elUni)
{

	//Sprawdzenie czy s¹ jakiekolwiek brzegowe granice
	bool isThereBC = false;
	for (int i = 0; i < 4; i++)
	{
		if (BC[i] == 1)
		{
			isThereBC = true;
			break;
		}
	}

	if (!isThereBC)
	{
		double** HBCForThisElement = new double* [4];

		for (int i = 0; i < 4; i++)
			HBCForThisElement[i] = new double[4];

		for (int n = 0; n < 4; n++)
		{
			for (int m = 0; m < 4; m++)
			{
				HBCForThisElement[n][m] = 0;
			}
		}

		HBC = HBCForThisElement;
	}

	double* G_X = initXGauss(elUni.nPkt);
	double* G_W = initWGauss(elUni.nPkt);

	std::vector<std::vector<double*>> edges;

	for (int i = 0; i < 4; i++)
	{
		std::vector<double*> edge;
		for (int j = 0; j < elUni.nPkt; j++)
		{
			double* pkt = new double[2];

			if (i % 2 == 0)
			{
				pkt[0] = G_X[elUni.nPkt - j - 1];
				pkt[1] = i != 0 ? 1 : -1;
			}
			else
			{
				pkt[0] = i != 1 ? -1 : 1;
				pkt[1] = G_X[elUni.nPkt - j - 1];
			}

			edge.push_back(pkt);
		}
		edges.push_back(edge);
	}


	std::cout << "Wspolrzedne punktow na krawedziach" << std::endl;

	for (int i = 0; i < edges.size(); i++)
	{
		for (int j = 0; j < elUni.nPkt; j++)
			std::cout << edges[i][j][0] << " " << edges[i][j][1] << std::endl;

		std::cout << std::endl;
	}


	std::vector<double**> matNPkt;


	for (int i = 0; i < edges.size(); i++)
	{
		double** N = new double* [elUni.nPkt];

		for (int i = 0; i < elUni.nPkt; i++)
			N[i] = new double[4];

		for (int j = 0; j < 2; j++)
		{

			N[j][0] = 0.25 * (1 - edges[i][j][0]) * (1 - edges[i][j][1]);
			N[j][1] = 0.25 * (1 + edges[i][j][0]) * (1 - edges[i][j][1]);
			N[j][2] = 0.25 * (1 + edges[i][j][0]) * (1 + edges[i][j][1]);
			N[j][3] = 0.25 * (1 - edges[i][j][0]) * (1 + edges[i][j][1]);
		}
		matNPkt.push_back(N);
	}

	std::cout << "Macierze funkcji N punktow" << std::endl;

	for (int i = 0; i < matNPkt.size(); i++)
	{
		for (int j = 0; j < elUni.nPkt; j++)
			std::cout << matNPkt[i][j][0] << " " << matNPkt[i][j][1] << " " << matNPkt[i][j][2] << " " << matNPkt[i][j][3] << std::endl;

		std::cout << std::endl;
	}

	std::vector<double> detJ;

	double L = sqrt(pow((x[0] - x[1]), 2) + pow((y[0] - y[1]), 2));
	detJ.push_back(L / 2);

	L = sqrt(pow((x[1] - x[2]), 2) + pow((y[1] - y[2]), 2));
	detJ.push_back(L / 2);

	L = sqrt(pow((x[2] - x[3]), 2) + pow((y[2] - y[3]), 2));
	detJ.push_back(L / 2);

	L = sqrt(pow((x[3] - x[0]), 2) + pow((y[3] - y[0]), 2));
	detJ.push_back(L / 2);

	std::vector<double**> HBCEdges;

	std::vector<int> edgesWithFlag;

	//Dol
	edgesWithFlag.push_back(BC[0] == 1 && BC[1] == 1 ? 1 : 0);
	edgesWithFlag.push_back(BC[1] == 1 && BC[2] == 1 ? 1 : 0);
	edgesWithFlag.push_back(BC[2] == 1 && BC[3] == 1 ? 1 : 0);
	edgesWithFlag.push_back(BC[3] == 1 && BC[0] == 1 ? 1 : 0);


	for (int p = 0; p < edges.size(); p++)
	{
		if (edgesWithFlag[p] == 0)
			continue;

		double** HBCEdge = new double* [4];

		for (int j = 0; j < 4; j++)
			HBCEdge[j] = new double[4];

		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
				HBCEdge[i][j] = 0;

			for (int j = 0; j < elUni.nPkt; j++)
			{
				for (int k = 0; k < 4; k++)
				{
					HBCEdge[i][k] += matNPkt[p][j][i] * matNPkt[p][j][k] * G_W[j] * 25.;
				}
			}

			for (int j = 0; j < 4; j++)
				HBCEdge[i][j] *= detJ[p];
		}

		HBCEdges.push_back(HBCEdge);


		std::cout << "Macierz HBC dla sciany " << p << std::endl;

		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
				std::cout << HBCEdge[i][j] << "\t";
			std::cout << std::endl;
		}

		std::cout << std::endl;

	}

	double** HBCForThisElement = new double* [4];

	for (int i = 0; i < 4; i++)
		HBCForThisElement[i] = new double[4];

	for (int n = 0; n < 4; n++)
	{
		for (int m = 0; m < 4; m++)
		{
			HBCForThisElement[n][m] = 0;
		}
	}

	for (int i = 0; i < HBCEdges.size(); i++)
	{
		for (int n = 0; n < 4; n++)
		{
			for (int m = 0; m < 4; m++)
			{
				HBCForThisElement[n][m] += HBCEdges[i][n][m];
			}
		}
	}

	std::cout << "Macierz HBC elementu" << std::endl;

	for (int n = 0; n < 4; n++)
	{
		for (int m = 0; m < 4; m++)
		{
			std::cout << HBCForThisElement[n][m] << " ";
		}
		std::cout << std::endl;
	}

	delete[] G_X;
	delete[] G_W;

	HBC = HBCForThisElement;
}

void Mesh::calcHForElements(ElementUniwersalny elUni)
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

void Mesh::calcHBCForElements(ElementUniwersalny elUni)
{
	for (int i = 0; i < globalData->ElementsNumber; i++)
	{
		std::cout << std::endl <<  "Element: " << i + 1 << std::endl;
		double* x = new double[4];
		double* y = new double[4];
		double* bc = new double[4];

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