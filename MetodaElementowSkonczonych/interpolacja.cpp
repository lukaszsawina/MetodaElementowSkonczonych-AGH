#include <iostream>
#include "interpolacja.h"

double Horner_PostacNaturalna(double* an, double x, int n)
{
	double wi = an[0];
	for (int i = 1; i <= n; i++)
		wi = wi * x + an[i];
	return wi;
}

double rekuHorner_PostacNewtona(double* an, double* xn, int x, int n)
{
	if (n == 0)
		return an[0];
	else
		return rekuHorner_PostacNewtona(an, xn, x, n - 1) * (x - xn[n - 1]) + an[n];
}

double newton_to_natural_recursive(double* bn, double* xn, int start, int end)
{
	if (start == end) {
		return bn[start];
	}
	else {
		double result = (newton_to_natural_recursive(bn, xn, start + 1, end) - newton_to_natural_recursive(bn, xn, start, end - 1)) / (xn[end] - xn[start]);
		bn[end] = result;
		return result;
	}
}

double* newton_to_natural(double* bn, double* xn, int n)
{
	for (int i = 1; i < n; i++) {
		newton_to_natural_recursive(bn, xn, 0, i);
	}
	return bn;
}

double lagrange(double* xn, double* yn, double x, int n)
{
	double output = 0;
	double ln;

	for (int i = 0; i < n; i++)
	{
		ln = 1;
		for (int j = 0; j < n; j++)
		{
			if (j == i)
				continue;
			ln *= (x - xn[j]) / (xn[i] - xn[j]);

		}
		output += yn[i] * ln;
	}
	return output;
}

double iloraz_roznicowy(double* xn, double* yn, int n)
{
	double wynik = 0;
	double mian;

	for (int i = 0; i < n; i++)
	{

		mian = 1;
		for (int j = 0; j < n; j++)
		{
			if (j == i)
				continue;
			mian *= (xn[i] - xn[j]);
		}
		wynik += (yn[i] / mian);
	}
	return wynik;
}