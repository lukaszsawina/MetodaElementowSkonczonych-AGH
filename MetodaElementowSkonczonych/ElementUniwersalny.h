#pragma once
#include <vector>

struct ElementUniwersalny {
	double** matdEta;
	double** matdKsi;
	std::vector<double**> matNPktForEdges;

	int nPkt;
	void init(int n);

	~ElementUniwersalny();
};
