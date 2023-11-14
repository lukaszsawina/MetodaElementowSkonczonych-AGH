#include <iostream>
#include <vector>

#include "calkowanie.h"
#include "MES.h"



//Zadanie domowe
double funkcja_testowa_1(double x) {
	// Definicja funkcji, któr¹ chcemy ca³kowaæ
	return 5 * pow(x,2) + 3 * x + 6; // Przyk³adowa funkcja kwadratowa
}

double funkcja_testowa_2(double x, double y) {
	// Definicja funkcji, któr¹ chcemy ca³kowaæ
	return 5 * pow(x,2) * pow(y,2) + 3 * x * y + 6; // Przyk³adowa funkcja kwadratowa
}


int main()
{
	//Mesh siatka;
	//siatka.readMeshFile("Test2_4_4_MixGrid.txt");
	//siatka.showGlobalData();
	//siatka.showNodes();
	//siatka.showElements();

	//ElementUniwersalny elementUniwersalny;
	//elementUniwersalny.init(2);

	//siatka.calcHForNodes(elementUniwersalny);

	//std::cout << "przestrzen 1d 2 punktowy schemat calkowania" << std::endl;
	//std::cout << Gauss1d(funkcja_testowa_1, -1, 1, 2) << std::endl;

	//std::cout << "przestrzen 1d 3 punktowy schemat calkowania" << std::endl;
	//std::cout << Gauss1d(funkcja_testowa_1, -1, 1, 3) << std::endl;

	//std::cout << "przestrzen 2d 2 punktowy schemat calkowania" << std::endl;
	//std::cout << Gauss2d(funkcja_testowa_2, -1, 1, 2) << std::endl;

	//std::cout << "przestrzen 2d 3 punktowy schemat calkowania" << std::endl;
	//std::cout << Gauss2d(funkcja_testowa_2, -1, 1, 3) << std::endl;


	ElementUniwersalny el;
	el.init(2);

	double x[4] = { 0, 0.025, 0.025, 0 };
	double y[4] = { 0, 0, 0.025, 0.025 };

	el.nPkt;

	double* G_X = initXGauss(el.nPkt);
	double* G_W = initWGauss(el.nPkt);

	std::vector<std::vector<double*>> edges;

	for (int i = 0; i < 4; i++)
	{
		std::vector<double*> edge;
		for (int j = 0; j < el.nPkt; j++)
		{
			double* pkt = new double[2];

			if (i % 2 == 0)
			{
				pkt[0] = G_X[el.nPkt - j - 1];
				pkt[1] = i != 0 ? 1 : -1;
			}
			else
			{
				pkt[0] = i != 1 ? -1 : 1;
				pkt[1] = G_X[el.nPkt - j - 1];
			}
			
			edge.push_back(pkt);
		}
		edges.push_back(edge);
	}


	std::cout << "Wspolrzedne punktow na krawedziach" << std::endl;

	for (int i = 0; i < edges.size(); i++)
	{
		for (int j = 0; j < el.nPkt; j++)
			std::cout << edges[i][j][0] << " " << edges[i][j][1] << std::endl;
	
		std::cout << std::endl;
	}


	std::vector<double**> matNPkt;
	

	for (int i = 0; i < edges.size(); i++)
	{
		double** N = new double* [el.nPkt];

		for (int i = 0; i < el.nPkt; i++)
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
		for (int j = 0; j < el.nPkt; j++)
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

	L = sqrt(pow((x[3] - x[1]), 2) + pow((y[3] - y[1]), 2));
	detJ.push_back(L / 2);

	std::vector<double**> HBCPkt;

	for (int p = 0; p < edges.size(); p++)
	{
		double** HBC = new double* [4];

		for (int j = 0; j < 4; j++)
			HBC[j] = new double[4];

		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
				HBC[i][j] = 0;

			for (int j = 0; j < el.nPkt; j++)
			{
				for (int k = 0; k < 4; k++)
				{
					HBC[i][k] += matNPkt[p][j][i] * matNPkt[p][j][k] * G_W[j] * 25.;
				}
			}

			for (int j = 0; j < 4; j++)
				HBC[i][j] *= detJ[p];
		}

		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
				std::cout << HBC[i][j] << "\t";
			std::cout << std::endl;
		}

		std::cout << std::endl;

	}

	


	delete[] G_X;
	delete[] G_W;


	return 0;
}

