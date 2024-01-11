#include "ElementUniwersalny.h"
#include "calkowanie.h"
#include "RunType.h"

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
	for (int i = 0, k = 0; i < n; i++)
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

	//Funkcje kszta³tu dla warunków brzegowych

	std::vector<std::vector<double*>> edges;

	for (int i = 0; i < 4; i++)
	{
		std::vector<double*> edge;
		for (int j = 0; j < nPkt; j++)
		{
			double* pkt = new double[2];

			if (i % 2 == 0)
			{
				pkt[0] = G_X[nPkt - j - 1];
				pkt[1] = i != 0 ? 1 : -1;
			}
			else
			{
				pkt[0] = i != 1 ? -1 : 1;
				pkt[1] = G_X[nPkt - j - 1];
			}

			edge.push_back(pkt);
		}
		edges.push_back(edge);
	}

	for (int i = 0; i < edges.size(); i++)
	{
		double** N = new double* [nPkt];

		for (int i = 0; i < nPkt; i++)
			N[i] = new double[4];

		for (int j = 0; j < nPkt; j++)
		{

			N[j][0] = 0.25 * (1 - edges[i][j][0]) * (1 - edges[i][j][1]);
			N[j][1] = 0.25 * (1 + edges[i][j][0]) * (1 - edges[i][j][1]);
			N[j][2] = 0.25 * (1 + edges[i][j][0]) * (1 + edges[i][j][1]);
			N[j][3] = 0.25 * (1 - edges[i][j][0]) * (1 + edges[i][j][1]);
		}
		matNPktForEdges.push_back(N);
	}

	std::vector<double*> pkt;

	for (int i = 0, k = 0; i < n; i++)
	{

		for (int j = 0; j < n; j++, k++)
		{
			double* Xpkt = new double[2];

			Xpkt[0] = G_X[nPkt - j - 1];
			Xpkt[1] = G_X[nPkt - i - 1];

			pkt.push_back(Xpkt);
		}
	}

	N = new double* [nPktCalk];

	for (int i = 0; i < nPktCalk; i++)
		N[i] = new double[4];

#ifdef DEBUG_ELEMENT_UNIWERSALNY
	std::cout << "--- Funkcje ksztaltow dla punktow ---" << std::endl;
#endif

	for (int j = 0; j < nPktCalk; j++)
	{
#ifdef DEBUG_ELEMENT_UNIWERSALNY
		std::cout << std::endl << "ksi: " << pkt[j][0] << " eta: " << pkt[j][1] << std::endl;
#endif

		N[j][0] = 0.25 * (1 - pkt[j][0]) * (1 - pkt[j][1]);
		N[j][1] = 0.25 * (1 + pkt[j][0]) * (1 - pkt[j][1]);
		N[j][2] = 0.25 * (1 + pkt[j][0]) * (1 + pkt[j][1]);
		N[j][3] = 0.25 * (1 - pkt[j][0]) * (1 + pkt[j][1]);


#ifdef DEBUG_ELEMENT_UNIWERSALNY
		std::cout << "N1: " << N[j][0];
		std::cout << " N2: " << N[j][1];
		std::cout << " N3: " << N[j][2];
		std::cout << " N4: " << N[j][3] << std::endl;
#endif

	}

	delete[] G_X;
}

ElementUniwersalny::~ElementUniwersalny()
{
	for (int i = 0; i <  pow(nPkt, 2); i++)
	{
		delete [] matdEta[i];
		delete [] matdKsi[i];
	}

	delete[] matdEta;
	delete[] matdKsi;

	std::vector<double**> matNPktForEdges;
	for (int i = 0; i < matNPktForEdges.size(); i++)
	{
		for (int j = 0; j < nPkt; j++)
			delete[] matNPktForEdges[i][j];

		delete[] matNPktForEdges[i];
	}
}