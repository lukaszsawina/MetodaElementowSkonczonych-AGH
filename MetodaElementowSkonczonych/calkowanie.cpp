#include "calkowanie.h"

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
	double nn = n;
	double* An = new double[nn];
	double* xn = new double[nn];

	switch (n)
	{
	case 2:
		xn[0] = -0.577350;
		xn[1] = 0.577350;
		An[0] = An[1] = 1;
		break;
	case 3:
		xn[0] = -0.774597;
		xn[1] = 0;
		xn[2] = 0.774597;
		An[0] = An[2] = (double)5 / 9;
		An[1] = (double)8 / 9;
		break;
	case 4:
		xn[0] = -0.861136;
		xn[1] = -0.339981;
		xn[2] = 0.339981;
		xn[3] = 0.861136;
		An[0] = An[3] = 0.347855;
		An[1] = An[2] = 0.652145;
		break;
	case 5:
		xn[0] = -0.906180;
		xn[1] = -0.538469;
		xn[2] = 0;
		xn[3] = 0.538469;
		xn[4] = 0.906180;
		An[0] = An[4] = 0.236927;
		An[1] = An[3] = 0.478629;
		An[2] = 0.568889;
		break;
	default:
		return 0;
	}

	for (int i = 0; i < nn; i++)
		output += An[i] * f(xn[i]);

	delete[] An;
	delete[] xn;
	return output;
}

double Gauss2d(double (*f)(double x, double y), double up, double down, int n)
{
	double output = 0.0;
	double nn = n;
	double* An = new double[nn];
	double* xn = new double[nn];

	switch (n)
	{
	case 2:
		xn[0] = -0.577350;
		xn[1] = 0.577350;
		An[0] = An[1] = 1;
		break;
	case 3:
		xn[0] = -0.774597;
		xn[1] = 0;
		xn[2] = 0.774597;
		An[0] = An[2] = (double)5 / 9;
		An[1] = (double)8 / 9;
		break;
	case 4:
		xn[0] = -0.861136;
		xn[1] = -0.339981;
		xn[2] = 0.339981;
		xn[3] = 0.861136;
		An[0] = An[3] = 0.347855;
		An[1] = An[2] = 0.652145;
		break;
	case 5:
		xn[0] = -0.906180;
		xn[1] = -0.538469;
		xn[2] = 0;
		xn[3] = 0.538469;
		xn[4] = 0.906180;
		An[0] = An[4] = 0.236927;
		An[1] = An[3] = 0.478629;
		An[2] = 0.568889;
		break;
	default:
		return 0;
	}

	for (int i = 0; i < nn; i++)
		for(int j = 0; j < nn; j++)
			output += An[i]*An[j] * f(xn[i], xn[j]);

	delete[] An;
	delete[] xn;
	return output;
}