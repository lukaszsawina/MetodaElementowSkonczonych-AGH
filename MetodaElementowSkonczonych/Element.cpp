#include <iostream>

#include "Element.h"
#include "calkowanie.h"
#include "RunType.h"

void Element::calcH(double* x, double* y, const ElementUniwersalny& elUni)
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

#ifdef DEBUG_ELEMENT_H
	std::cout << "---Macierze jakobiego punktow---" << std::endl;
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
#endif


	//Przemnożenie przez 1/detJ
	for (int i = 0; i < macierzeJakobiegoPunktow.size(); i++)
	{
		for (int j = 0; j < 2; j++)
			for (int k = 0; k < 2; k++)
				macierzeJakobiegoPunktow[i][j][k] *= (1 / detJPunktow[i]);
	}

#ifdef DEBUG_ELEMENT_H
	std::cout << std::endl << "---Macierze jakobiego punktów przemnożona przez 1/detJ---" << std::endl;
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
#endif

	double** matdx;
	double** matdy;

	matdx = new double* [pow(elUni.nPkt, 2)];
	matdy = new double* [pow(elUni.nPkt, 2)];

	for (int i = 0; i < pow(elUni.nPkt, 2); i++)
	{
		matdx[i] = new double[4];
		matdy[i] = new double[4];
	}

	//Liczenie wartości macierzy dx i dy
	for (int i = 0; i < pow(elUni.nPkt, 2); i++)
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

	//Zwolnienei pamięci z macierzyJakobiegoPunktów
	for (int i = 0; i < macierzeJakobiegoPunktow.size(); i++)
	{
		for (int j = 0; j < 2; j++)
			delete[] macierzeJakobiegoPunktow[i][j];

		delete[] macierzeJakobiegoPunktow[i];
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
			Hpkt[i][0] = (matdx[p][i] * matdx[p][0] + matdy[p][i] * matdy[p][0]) * detJPunktow[p] * globalData->Conductivity;
			Hpkt[i][1] = (matdx[p][i] * matdx[p][1] + matdy[p][i] * matdy[p][1]) * detJPunktow[p] * globalData->Conductivity;
			Hpkt[i][2] = (matdx[p][i] * matdx[p][2] + matdy[p][i] * matdy[p][2]) * detJPunktow[p] * globalData->Conductivity;
			Hpkt[i][3] = (matdx[p][i] * matdx[p][3] + matdy[p][i] * matdy[p][3]) * detJPunktow[p] * globalData->Conductivity;
		}

		macierzeHPunktow.push_back(Hpkt);
	}


	//Zwalnianie pamięci z matdx matdy
	for (int i = 0; i < 4; i++)
	{
		delete[] matdx[i];
		delete[] matdy[i];
	}

	delete[] matdx;
	delete[] matdy;

#ifdef DEBUG_ELEMENT_H
	std::cout << std::endl;
	std::cout << "---Macierze H punktow---" << std::endl;

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
#endif

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

	//Zwolnienie pamięci dla macierzy H punktów
	delete[] G_W;
	for (int i = 0; i < macierzeHPunktow.size(); i++)
	{
		for (int j = 0; j < 4; j++)
			delete[] macierzeHPunktow[i][j];

		delete[] macierzeHPunktow[i];
	}
#ifdef DEBUG_ELEMENT_H
	std::cout << "Macierz H" << std::endl;
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
			std::cout << outputH[i][j] << "\t";
		std::cout << std::endl;
	}
#endif

	H = outputH;
}

void Element::calcHBC(double* x, double* y, int* BC, const ElementUniwersalny& elUni)
{

	//Sprawdzenie czy są jakiekolwiek brzegowe granice
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


	for (int p = 0; p < 4; p++)
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
					HBCEdge[i][k] += elUni.matNPktForEdges[p][j][i] * elUni.matNPktForEdges[p][j][k] * G_W[j] * globalData->Alfa;
				}
			}

			for (int j = 0; j < 4; j++)
				HBCEdge[i][j] *= detJ[p];
		}

		HBCEdges.push_back(HBCEdge);

#ifdef DEBUG_ELEMENT_HBC
		std::cout << "Macierz HBC dla sciany " << p << std::endl;

		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
				std::cout << HBCEdge[i][j] << "\t";
			std::cout << std::endl;
		}

		std::cout << std::endl;
#endif
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

	//Zwolnienie pamieci z HBCEdges

	for (int i = 0; i < HBCEdges.size(); i++)
	{
		for (int j = 0; j < 4; j++)
		{
			delete[] HBCEdges[i][j];
		}
		delete[] HBCEdges[i];
	}

#ifdef DEBUG_ELEMENT_HBC
	std::cout << "Macierz HBC elementu" << std::endl;

	for (int n = 0; n < 4; n++)
	{
		for (int m = 0; m < 4; m++)
		{
			std::cout << HBCForThisElement[n][m] << " ";
		}
		std::cout << std::endl;
	}
#endif

	delete[] G_X;
	delete[] G_W;

	HBC = HBCForThisElement;
}

void Element::calcVectorP(double* x, double* y, int* BC, const ElementUniwersalny& elUni)
{
	double* G_W = initWGauss(elUni.nPkt);

	std::vector<int> edgesWithFlag;

	//Dol
	edgesWithFlag.push_back(BC[0] == 1 && BC[1] == 1 ? 1 : 0);
	edgesWithFlag.push_back(BC[1] == 1 && BC[2] == 1 ? 1 : 0);
	edgesWithFlag.push_back(BC[2] == 1 && BC[3] == 1 ? 1 : 0);
	edgesWithFlag.push_back(BC[3] == 1 && BC[0] == 1 ? 1 : 0);

	std::vector<double> detJ;

	double L = sqrt(pow((x[0] - x[1]), 2) + pow((y[0] - y[1]), 2));
	detJ.push_back(L / 2);

	L = sqrt(pow((x[1] - x[2]), 2) + pow((y[1] - y[2]), 2));
	detJ.push_back(L / 2);

	L = sqrt(pow((x[2] - x[3]), 2) + pow((y[2] - y[3]), 2));
	detJ.push_back(L / 2);

	L = sqrt(pow((x[3] - x[0]), 2) + pow((y[3] - y[0]), 2));
	detJ.push_back(L / 2);

	std::vector<double*> vectorPForEdges;

	for (int e = 0; e < 4; e++)
	{
		double* vectorPForEdge = new double[4];

		for (int i = 0; i < 4; i++)
			vectorPForEdge[i] = 0;

		if (edgesWithFlag[e] == 0)
		{
			vectorPForEdges.push_back(vectorPForEdge);
			continue;
		}

		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < elUni.nPkt; j++)
			{
				vectorPForEdge[i] += elUni.matNPktForEdges[e][j][i] * globalData->Tot * G_W[j];
			}
		}

		for (int j = 0; j < 4; j++)
			vectorPForEdge[j] = vectorPForEdge[j] * detJ[e] * globalData->Alfa;

		vectorPForEdges.push_back(vectorPForEdge);
	}

#ifdef DEBUG_ELEMENT_VP
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
			std::cout << vectorPForEdges[i][j] << "\t";
		std::cout << std::endl;
	}
#endif

	double* wynik = new double[4];

	for (int i = 0; i < 4; i++)
		wynik[i] = 0;

	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
			wynik[i] += vectorPForEdges[j][i];
	}

	//Zwolnienie pamięci z vectorPForEdges
	for (int i = 0; i < vectorPForEdges.size(); i++)
	{
		delete[] vectorPForEdges[i];
	}

#ifdef DEBUG_ELEMENT_VP
	for (int i = 0; i < 4; i++)
	{
		std::cout << wynik[i] << "\t";
	}
	std::cout << std::endl;
#endif

	VectorP = wynik;
}

void Element::calcC(double* x, double* y, const ElementUniwersalny& elUni)
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

	//Przemnożenie przez 1/detJ
	for (int i = 0; i < macierzeJakobiegoPunktow.size(); i++)
	{
		for (int j = 0; j < 2; j++)
			for (int k = 0; k < 2; k++)
				macierzeJakobiegoPunktow[i][j][k] *= (1 / detJPunktow[i]);
	}

	double** outputC = new double* [4];
	for (int i = 0; i < 4; i++)
		outputC[i] = new double[4];

	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
			outputC[i][j] = 0;
	}

	std::vector<double**> macierzeCPunktow;

	for (int p = 0; p < pow(elUni.nPkt, 2); p++)
	{
		double** Cpkt = new double* [4];

		for (int i = 0; i < 4; i++)
			Cpkt[i] = new double[4];

		//Liczenie macierzy C dla pc[i]
		for (int i = 0; i < 4; i++)
		{
			Cpkt[i][0] = (elUni.N[p][i] * elUni.N[p][0]) * detJPunktow[p] * globalData->Density * globalData->SpecificHeat;
			Cpkt[i][1] = (elUni.N[p][i] * elUni.N[p][1]) * detJPunktow[p] * globalData->Density * globalData->SpecificHeat;
			Cpkt[i][2] = (elUni.N[p][i] * elUni.N[p][2]) * detJPunktow[p] * globalData->Density * globalData->SpecificHeat;
			Cpkt[i][3] = (elUni.N[p][i] * elUni.N[p][3]) * detJPunktow[p] * globalData->Density * globalData->SpecificHeat;
		}

		macierzeCPunktow.push_back(Cpkt);
	}

	double* G_W = initWGauss(elUni.nPkt);

	//Liczenie macierzy C dla elementu
	for (int i = 0, p = 0; i < elUni.nPkt; i++)
	{
		for (int j = 0; j < elUni.nPkt; j++, p++)
		{
			for (int n = 0; n < 4; n++)
				for (int m = 0; m < 4; m++)
				{
					macierzeCPunktow[p][n][m] *= (G_W[i] * G_W[j]);
					outputC[n][m] += macierzeCPunktow[p][n][m];
				}
		}
	}

#ifdef DEBUG_ELEMENT_C
	std::cout << "Macierz C" << std::endl;
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
			std::cout << outputC[i][j] << "\t";
		std::cout << std::endl;
	}
#endif

	C = outputC;
}