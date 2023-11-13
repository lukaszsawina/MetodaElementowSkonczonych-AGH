#include "calkowanie.h"


double* initXGauss(int n)
{
	double* G_X = new double[n];

	switch (n)
	{
	case 2:
		G_X[0] = -sqrt(1./3.);
		G_X[1] = sqrt(1./3.);
		break;
	case 3:
		G_X[0] = -sqrt(3./5.);
		G_X[1] = 0;
		G_X[2] = sqrt(3./5.);
		break;
	case 4:
		G_X[0] = -sqrt((3. / 7.) + (2. / 7.) * sqrt(6. / 5.));
		G_X[1] = -sqrt((3. / 7.) - (2. / 7.) * sqrt(6. / 5.));
		G_X[2] = sqrt((3. / 7.) - (2. / 7.) * sqrt(6. / 5.));
		G_X[3] = sqrt((3. / 7.) + (2. / 7.) * sqrt(6. / 5.));
		break;
	case 5:
		G_X[0] = -(1. / 3.) * sqrt(5. + 2. * (10. / 7.));
		G_X[1] = -(1. / 3.) * sqrt(5. - 2. * (10. / 7.));
		G_X[2] = 0;
		G_X[3] = (1. / 3.) * sqrt(5. - 2. * (10. / 7.));
		G_X[4] = (1. / 3.) * sqrt(5. + 2. * (10. / 7.));
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
		G_W[0] = G_W[2] = (5. / 9.);
		G_W[1] = (8. / 9.);
		break;
	case 4:
		G_W[0] = G_W[3] = (18. - sqrt(30.))/36.;
		G_W[1] = G_W[2] = (18. + sqrt(30.)) / 36.;
		break;
	case 5:
		G_W[0] = G_W[4] = (322.-13.*sqrt(70.))/900.;
		G_W[1] = G_W[3] = (322. + 13. * sqrt(70.)) / 900.;
		G_W[2] = (128./225.);
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