#pragma once
#include <vector>

struct ElementUniwersalny {
	double** matdEta;
	double** matdKsi;
	std::vector<double**> matNPktForEdges;
	double** N;

	int nPkt;
	void init(int n);

	~ElementUniwersalny();
};
