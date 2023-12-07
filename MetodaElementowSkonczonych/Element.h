#pragma once
#include "ElementUniwersalny.h"

struct Element {
	int ID;
	int ID_wezlow[4];
	double** H = nullptr;
	double** HBC = nullptr;
	double* VectorP = nullptr;
	double** C = nullptr;

	void calcH(double* x, double* y, const ElementUniwersalny& elUni);
	void calcHBC(double* x, double* y, int* BC, const ElementUniwersalny& elUni);
	void calcVectorP(double* x, double* y, int* BC, const ElementUniwersalny& elUni);
	void calcC(double* x, double* y, const ElementUniwersalny& elUni);
};