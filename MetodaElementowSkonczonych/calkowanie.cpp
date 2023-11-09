#include "calkowanie.h"


double* initXGauss(int n)
{
	double* G_X = new double[n];

	switch (n)
	{
	case 2:
		G_X[0] = -0.577350;
		G_X[1] = 0.577350;
		break;
	case 3:
		G_X[0] = -0.774597;
		G_X[1] = 0;
		G_X[2] = 0.774597;
		break;
	case 4:
		G_X[0] = -0.861136;
		G_X[1] = -0.339981;
		G_X[2] = 0.339981;
		G_X[3] = 0.861136;
		break;
	case 5:
		G_X[0] = -0.906180;
		G_X[1] = -0.538469;
		G_X[2] = 0;
		G_X[3] = 0.538469;
		G_X[4] = 0.906180;
		break;
	default:
		return 0;
	}

	return G_X;
}

double* initWGauss(int n)
{
	double* G_W = new double[n];

	switch (n)
	{
	case 2:
		G_W[0] = G_W[1] = 1;
		break;
	case 3:
		G_W[0] = G_W[2] = (double)5 / 9;
		G_W[1] = (double)8 / 9;
		break;
	case 4:
		G_W[0] = G_W[3] = 0.347855;
		G_W[1] = G_W[2] = 0.652145;
		break;
	case 5:
		G_W[0] = G_W[4] = 0.236927;
		G_W[1] = G_W[3] = 0.478629;
		G_W[2] = 0.568889;
		break;
	default:
		return 0;
	}

	return G_W;
}


double Calkowanie_metoda_prostokatow(double (*f)(double x), double a, double b, int n)
{
	double output = 0;
	double offset = fabs(b - a) / n;
	double dx = offset / 2;
	double j = a;


	for (int i = 0; i < n; i++)
	{
		output += offset * f(j);
		j += offset;
	}

	return output;
}

double Gauss1d(double (*f)(double x), double up, double down, int n)
{
	double output = 0.0;
	double* G_W = initWGauss(n);
	double* G_X = initXGauss(n);

	for (int i = 0; i < n; i++)
		output += G_W[i] * f(G_X[i]);

	delete[] G_W;
	delete[] G_X;
	return output;
}

double Gauss2d(double (*f)(double x, double y), double up, double down, int n)
{
	double output = 0.0;
	double* G_W = initWGauss(n);
	double* G_X = initXGauss(n);

	for (int i = 0; i < n; i++)
		for(int j = 0; j < n; j++)
			output += G_W[i]* G_W[j] * f(G_X[i], G_X[j]);

	delete[] G_W;
	delete[] G_X;
	return output;
}