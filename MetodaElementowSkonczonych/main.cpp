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
    //GlobalData* data = readMesh("Test1_4_4.txt");
	//data->showGlobalData();

	//std::cout << "przestrzen 1d 2 punktowy schemat calkowania" << std::endl;
	//std::cout << Gauss1d(funkcja_testowa_1, -1, 1, 2) << std::endl;

	//std::cout << "przestrzen 1d 3 punktowy schemat calkowania" << std::endl;
	//std::cout << Gauss1d(funkcja_testowa_1, -1, 1, 3) << std::endl;

	//std::cout << "przestrzen 2d 2 punktowy schemat calkowania" << std::endl;
	//std::cout << Gauss2d(funkcja_testowa_2, -1, 1, 2) << std::endl;

	//std::cout << "przestrzen 2d 3 punktowy schemat calkowania" << std::endl;
	//std::cout << Gauss2d(funkcja_testowa_2, -1, 1, 3) << std::endl;

	ElementUniwersalny elementUniwersalny;

	elementUniwersalny.init(2);

	double x[4] = { 0, 0.025, 0.025, 0 };
	double y[4] = { 0, 0, 0.025, 0.025 };

	std::vector<double**> macierzeJakobiegoPunktow;

	std::vector<double> detJPunktow;
	
	for (int i = 0; i < pow(elementUniwersalny.nPkt,2); i++)
	{
		double** macierzJakobiego = new double* [2];

		for (int i = 0; i < 2; i++)
		{
			macierzJakobiego[i] = new double[2];
		}

		//Liczenie macierzy jakobiego dla pc[i]
		macierzJakobiego[0][0] = elementUniwersalny.matdEta[i][0] * y[0] + elementUniwersalny.matdEta[i][1] * y[1] + elementUniwersalny.matdEta[i][2] * y[2] + elementUniwersalny.matdEta[i][3] * y[3];
		macierzJakobiego[0][1] = -1 * (elementUniwersalny.matdKsi[i][0] * y[0] + elementUniwersalny.matdKsi[i][1] * y[1] + elementUniwersalny.matdKsi[i][2] * y[2] + elementUniwersalny.matdKsi[i][3] * y[3]);
		macierzJakobiego[1][0] = -1 * (elementUniwersalny.matdEta[i][0] * x[0] + elementUniwersalny.matdEta[i][1] * x[1] + elementUniwersalny.matdEta[i][2] * x[2] + elementUniwersalny.matdEta[i][3] * x[3]);
		macierzJakobiego[1][1] = elementUniwersalny.matdKsi[i][0] * x[0] + elementUniwersalny.matdKsi[i][1] * x[1] + elementUniwersalny.matdKsi[i][2] * x[2] + elementUniwersalny.matdKsi[i][3] * x[3];

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
			for(int k = 0; k < 2; k++)
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

	for (int i = 0; i < 4; i++)
	{
		matdx[i][0] = macierzeJakobiegoPunktow[i][0][0] * elementUniwersalny.matdKsi[i][0] + macierzeJakobiegoPunktow[i][0][1] * elementUniwersalny.matdEta[i][0];
		matdx[i][1] = macierzeJakobiegoPunktow[i][0][0] * elementUniwersalny.matdKsi[i][1] + macierzeJakobiegoPunktow[i][0][1] * elementUniwersalny.matdEta[i][1];
		matdx[i][2] = macierzeJakobiegoPunktow[i][0][0] * elementUniwersalny.matdKsi[i][2] + macierzeJakobiegoPunktow[i][0][1] * elementUniwersalny.matdEta[i][2];
		matdx[i][3] = macierzeJakobiegoPunktow[i][0][0] * elementUniwersalny.matdKsi[i][3] + macierzeJakobiegoPunktow[i][0][1] * elementUniwersalny.matdEta[i][3];

		matdy[i][0] = macierzeJakobiegoPunktow[i][1][0] * elementUniwersalny.matdKsi[i][0] + macierzeJakobiegoPunktow[i][1][1] * elementUniwersalny.matdEta[i][0];
		matdy[i][1] = macierzeJakobiegoPunktow[i][1][0] * elementUniwersalny.matdKsi[i][1] + macierzeJakobiegoPunktow[i][1][1] * elementUniwersalny.matdEta[i][1];
		matdy[i][2] = macierzeJakobiegoPunktow[i][1][0] * elementUniwersalny.matdKsi[i][2] + macierzeJakobiegoPunktow[i][1][1] * elementUniwersalny.matdEta[i][2];
		matdy[i][3] = macierzeJakobiegoPunktow[i][1][0] * elementUniwersalny.matdKsi[i][3] + macierzeJakobiegoPunktow[i][1][1] * elementUniwersalny.matdEta[i][3];
	}

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

	for (int p = 0; p < pow(elementUniwersalny.nPkt, 2); p++)
	{
		double** Hpkt = new double* [4];

		for (int i = 0; i < 4; i++)
			Hpkt[i] = new double[4];

		for (int i = 0; i < 4; i++)
		{
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

	
	int nn = elementUniwersalny.nPkt;
	double* An = new double[nn];

	switch (nn)
	{
	case 2:
		An[0] = An[1] = 1;
		break;
	case 3:
		An[0] = An[2] = (double)5 / 9;
		An[1] = (double)8 / 9;
		break;
	case 4:
		An[0] = An[3] = 0.347855;
		An[1] = An[2] = 0.652145;
		break;
	case 5:
		An[0] = An[4] = 0.236927;
		An[1] = An[3] = 0.478629;
		An[2] = 0.568889;
		break;
	}

	for (int i = 0, p = 0; i < elementUniwersalny.nPkt; i++)
	{
		for (int j = 0; j < elementUniwersalny.nPkt; j++, p++)
		{
			for (int n = 0; n < 4; n++)
				for (int m = 0; m < 4; m++)
				{
					macierzeHPunktow[p][n][m] *= (An[i] * An[j]);
					H[n][m] += macierzeHPunktow[p][n][m];
				}
		}
	}

	std::cout << "Macierz H" << std::endl;
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
			std::cout << H[i][j] << "\t";
		std::cout << std::endl;
	}



	return 0;
}

